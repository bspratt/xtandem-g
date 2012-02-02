/*
 * Copyright (c) 2012 Insilicos LLC
 *
 * A CUDA-based reimplementation of work which is
 * Copyright (c) 2003-2006 Fred Hutchinson Cancer Research Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <iostream>
#include <sstream>

#include "mscore_kgpu_thrust.h"
#include <cuda.h>
#include <device_functions.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/transform_reduce.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>
#include <thrust/unique.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>

typedef thrust::device_vector<float> fvec;
typedef thrust::device_vector<int> ivec;
typedef ivec::iterator iveciter;
typedef fvec::iterator fveciter;

void vmiTypeGPU::init(int s) {
    fI = new thrust::device_vector<float>(s);
    iM = new thrust::device_vector<int>(s);
}
void vmiTypeGPU::kill() {
    delete fI;
    fI = NULL;
    delete iM;
    iM = NULL;
}

// helpers for hiding device memory implementation from normal C++
thrust::device_vector<float> *mscore_kgpu_thrust_fvec_alloc(int size) {
    return new thrust::device_vector<float>(size);
}
void mscore_kgpu_thrust_fvec_clear(thrust::device_vector<float> *vec) {
    thrust::fill(vec->begin(), vec->end(), 0);
}
void mscore_kgpu_thrust_fvec_kill(thrust::device_vector<float> *vec) {
    delete vec;
}


// define transformation f(x) -> x^2
struct square
{
    __host__ __device__
        float operator()(float x)
    {
        return x * x;
    }
};

struct divideby
{
    const float m_factor;
    __host__ __device__
        divideby(double divisor) :  m_factor ((float)(1.0/divisor)) {}
    __host__ __device__
        float operator()(float m)
    {
        return (m*m_factor);
    }
};

struct multiplyby
{
    const float m_factor;
    __host__ __device__
        multiplyby(double factor) : m_factor(factor) {}
    __host__ __device__
        float operator()(float m)
    {
        return (m*m_factor);
    }
};

struct incrementby
{
    const int m_incr;
    __host__ __device__
        incrementby(int incr) :  m_incr (incr) {}
    __host__ __device__
        int operator()(int m)
    {
        return (m+m_incr);
    }
};


/*
* m/z binning
*/
struct imass_functor
{
    const float m_IsotopeCorrection;
    const int m_binOffset;
    __host__ __device__
        imass_functor(double IsotopeCorrection, int binOffset) : 
    m_IsotopeCorrection((float)(1.0/IsotopeCorrection)),m_binOffset(binOffset) {}
    template <typename Tuple>
    __host__ __device__
        void operator()(Tuple m)
    {
        thrust::get<1>(m) =  ((int)((thrust::get<0>(m)*m_IsotopeCorrection) + 0.5f))-m_binOffset;
    }
};


/*
* mixrange
* fI[i] = tempLookup[i] - (sums[i+101]-sums[i])/101
* <0> = <1> - (<2>-<3>)/101
*/
struct mixrange_functor
{
    template <typename Tuple>
    __host__ __device__
        void operator()(Tuple m) 
    {
        thrust::get<0>(m) = thrust::get<1>(m) -  ((thrust::get<2>(m)-thrust::get<3>(m))/101.0f);
    }
};

/*
* resolve mass binning collisions
*/
struct imass_collision_functor {
    template <typename Tuple>
    __host__ __device__
        void operator()(Tuple m) // mass a, mass b, intensity a, intensity b
    {
        if (thrust::get<1>(m) == thrust::get<0>(m)) { // bin collision
            if (thrust::get<2>(m) < thrust::get<3>(m)) { // take max intensity
                thrust::get<2>(m) = thrust::get<3>(m);
            } else {
                thrust::get<3>(m) = thrust::get<2>(m);
            }
        }
    }
};

struct nonpositive_functor
{
    template <typename Tuple>
    __host__ __device__
        bool operator()(Tuple m) 
    {
        return (thrust::get<1>(m) <= 0);
    }
};

// define transformation f(x) -> sqrt(x)
struct squareroot
{
    __host__ __device__
        float operator()(float x)
    {
        return sqrt(x);
    }
};


void mscore_kgpu_thrust_score(thrust::host_vector<float> const &host_fI, // intensities (input)
    thrust::host_vector<float> const &host_fM, // m/z's (input)
    double dIsotopeCorrection, // for m/z binning (input)
    int iWindowCount, // (input)
    int endMassMax, // max acceptable binned m/z (input)
    int &m_maxEnd,  // max encountered binned m/z (output)
    vmiTypeGPU &binned)   // binned intensities and m/z's (output)
{
    if (!host_fM.size()) {
        return;
    }
    thrust::device_vector<int> &iM = *binned.iM;
    thrust::device_vector<float> &fI = *binned.fI;

    // figure the shift needed for lowest binned mass and half smoothing window
    int startMass = (int)((host_fM[0]*dIsotopeCorrection) + 0.5f);
    const int binOffset = startMass-50;
    endMassMax -= binOffset; // shift to local coords
    startMass -= binOffset;

    // bin the m/z's
    iM.resize(host_fM.size());
    {
    fvec fM(host_fM); // copy to device memory
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(fM.begin(), iM.begin())), 
                     thrust::make_zip_iterator(thrust::make_tuple(fM.end(), iM.end())), 
                     imass_functor(dIsotopeCorrection,binOffset));
    }

    // Screen peeks on upper end.  TODO faster way?
    iveciter itM = iM.begin();
    iveciter itEnd = iM.end();
    int endMass = itEnd[-1];
    while (itM != itEnd && endMass >= endMassMax) {
        itEnd--;
        endMass = itEnd[-1];
    }

    if (itM == itEnd)  // No peaks left.
    {
        return;
    }

    // if any bin has more than one occupant, choose the one with higher intensity
    fI = host_fI; // copy to device memory
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(iM.begin(), iM.begin()+1, fI.begin(), fI.begin()+1)), 
                     thrust::make_zip_iterator(thrust::make_tuple(iM.end()-1, iM.end(), fI.end()-1, fI.end())), 
                     imass_collision_functor());

    // take sqrt of intensities
    thrust::transform(fI.begin(),fI.end(),fI.begin(),squareroot());

    fveciter itI = fI.begin();

    // for every pair fI,iM set templookup[iM]=fI
    fvec tempLookup(endMass+50); // +50 for smoothing window
    thrust::scatter(fI.begin(),fI.end(),iM.begin(),tempLookup.begin());
    if ((endMass+binOffset)+50 > m_maxEnd)
        m_maxEnd = (endMass+binOffset)+50;

    float fMaxI = *thrust::max_element(fI.begin(),fI.end());
    float fMinCutoff = (float) (0.05 * fMaxI);

    int range = thrust::min(endMassMax, iWindowCount + endMass)-startMass;
    if (range > 3000)
        iWindowCount=10;
    else if (range > 2500)
        iWindowCount=9;
    else if (range > 2000)
        iWindowCount=8;
    else if (range > 1500)
        iWindowCount=7;
    else if (range > 1000)
        iWindowCount=6;
    else
        iWindowCount=5;

    int iWindowSize = range / iWindowCount;

    /*
     * Process input spectrum for dot product - split windows
     */
    for (int i = 0; i < iWindowCount; i++) {

        int iStart = startMass + i*iWindowSize;

        /*
         * Get maximum intensity within window
         */
        float fMaxWindowI = *thrust::max_element(tempLookup.begin()+iStart,tempLookup.begin()+(iStart+iWindowSize));

        if (fMaxWindowI > 0.0 && fMaxWindowI > fMinCutoff) {
            /*
             * Normalize within window
             */
            thrust::transform(tempLookup.begin()+iStart,tempLookup.begin()+(iStart + iWindowSize),
                tempLookup.begin()+iStart,multiplyby(fMaxI / fMaxWindowI));
        }
    }

 /*
    * Reduce intensity and make unit vector by dividing
    * every point by sqrt(sum(x^2))
    */
    double dSpectrumArea = sqrt(thrust::transform_reduce(tempLookup.begin()+startMass, tempLookup.begin()+(endMass+1), square(), 0.0f, thrust::plus<float>()));
    thrust::transform(tempLookup.begin()+startMass,tempLookup.begin()+(endMass+1),tempLookup.begin()+startMass,divideby(dSpectrumArea));

    /*
     * Perform mix-range modification to input spectrum
     */
    thrust::device_vector<float> sums(tempLookup.size()+101);
    thrust::copy(tempLookup.begin(),tempLookup.end(),sums.begin()+50);
    thrust::exclusive_scan(sums.begin()+50, sums.end(), sums.begin()+50);
    // now sums[i+101]-sums[i] = sum(tempLookup[i-50:i+50])
    fI.resize(tempLookup.size());
    // fI[i] = tempLookup[i] - (sums[i+101]-sums[i])/101
    thrust::for_each(thrust::make_zip_iterator(
                     thrust::make_tuple( fI.begin(), tempLookup.begin(), sums.begin()+101, sums.begin())),
                     thrust::make_zip_iterator(
                     thrust::make_tuple( fI.end(), tempLookup.end(), sums.end(),sums.end()-101)), 
                     mixrange_functor());

    // now remove any non-positive results
    iM.resize(fI.size());
    thrust::sequence(iM.begin(),iM.end(),binOffset,1); // shift back to host coords
     thrust::zip_iterator<thrust::tuple<iveciter,fveciter>> new_end = 
       thrust::remove_if(thrust::make_zip_iterator(
        thrust::make_tuple(iM.begin(),fI.begin())),
        thrust::make_zip_iterator(
        thrust::make_tuple(iM.end(),fI.end())),
        nonpositive_functor());
    // trim to new length
    iM.erase(thrust::get<0>(new_end.get_iterator_tuple()),iM.end());
    fI.erase(thrust::get<1>(new_end.get_iterator_tuple()),fI.end());
}


/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses).
 */

struct dot_functor
{
    __host__ __device__
        dot_functor() {}
    template <typename Tuple>
    __host__ __device__
        void operator()(Tuple m) // 0=sequence ion, 1=spectrum ion, 2=intensity, 3=used, 4=dscore, 5=lcount
    {
        if (thrust::get<0>(m) == thrust::get<1>(m)) { // ion found in spectrum and sequence
            if (thrust::get<3>(m) < thrust::get<2>(m)) { // more intense than previous m_used, if any
                thrust::get<4>(m) += thrust::get<2>(m) - thrust::get<3>(m); // dscore
                thrust::get<3>(m) = thrust::get<2>(m); // m_used
                thrust::get<5>(m)++; // lcount
            }
        }
    }
};
struct dot_neighbor_functor
{
    __host__ __device__
        dot_neighbor_functor() {}
    template <typename Tuple>
    __host__ __device__
        void operator()(Tuple m) // 0=sequence ion, 1=spectrum ion, 2=spectrum ion intensity, 3=spectrum ion used, 4=dscore
    {
        if ((thrust::get<0>(m)) == thrust::get<1>(m)) { // neighbor is right next door
            float fTmp = 0.5f*thrust::get<2>(m); // use 1/2 intensity
            if (thrust::get<3>(m) < fTmp) { // more intense than previous m_used, if any
                thrust::get<4>(m) += fTmp - thrust::get<3>(m); // dscore
                thrust::get<3>(m) = fTmp; // m_used
            }
        }
    }
};
struct add_B_to_A_functor
{
    template <typename Tuple>
    __host__ __device__
        void operator()(Tuple m) 
    {
        thrust::get<0>(m) += thrust::get<1>(m);
    }
};

// helpful macro for copying to c++ space for debugger viewing
#define STDVECT(T,foo,seq) std::vector<T> foo; for (int n=0;n<seq.size();foo.push_back(seq[n++]));

float mscore_kgpu_thrust_dot(unsigned long &lCount,const vmiTypeGPU &_spectrum,unsigned long const *_plSeq,
                    thrust::device_vector<float> &_miUsed) {
    thrust::host_vector<int> hseq;  // TODO get this into device memory sooner (if you do, see note "restore" below
    for (unsigned long const *p = _plSeq;*p;) {
        hseq.push_back(*p++);
    }

    thrust::device_vector<int> seq(hseq);

    // for each peak mz in current sequence, find the first peak in the spectrum
    // of greater or equal mz, store map in seqindex
    thrust::device_vector<int> seqindex(seq.size()); 
    thrust::lower_bound(_spectrum.iM->begin(),_spectrum.iM->end(),seq.begin(),seq.end(),seqindex.begin());

    // now do dot on peak mz values that appear in seq and spectrum
    thrust::device_vector<float> dscores(seq.size());
    thrust::device_vector<int> lcounts(seq.size());
    thrust::for_each(
        thrust::make_zip_iterator( 
            thrust::make_tuple(
                seq.begin(),
                thrust::make_permutation_iterator(_spectrum.iM->begin(), seqindex.begin() ), 
                thrust::make_permutation_iterator(_spectrum.fI->begin(), seqindex.begin() ), 
                thrust::make_permutation_iterator(_miUsed.begin(), seq.begin() ),
                dscores.begin(),
                lcounts.begin())), 
        thrust::make_zip_iterator(
            thrust::make_tuple(
                seq.end(),
                thrust::make_permutation_iterator(_spectrum.iM->begin(), seqindex.end()), 
                thrust::make_permutation_iterator(_spectrum.fI->begin(), seqindex.end() ), 
                thrust::make_permutation_iterator(_miUsed.begin(), seq.end() ),
                dscores.end(),
                lcounts.end())), 
         dot_functor());

    lCount = thrust::reduce(lcounts.begin(),lcounts.end()); // sum of all lcounts

    // now check for left and right neighbors to contribute 50%
    // NOTE: seems whacky that we would do the left-right thing even when 
    // the actual ion did not match, but that's how k-score does it - bpratt
    thrust::device_vector<int> neighbors(seq);
    // look left
    thrust::transform(neighbors.begin(),neighbors.end(),neighbors.begin(),incrementby(-1));
    thrust::transform(seqindex.begin(),seqindex.end(),seqindex.begin(),incrementby(-1));
    thrust::for_each(
        thrust::make_zip_iterator( 
            thrust::make_tuple(
                neighbors.begin(), 
                thrust::make_permutation_iterator(_spectrum.iM->begin(), seqindex.begin()), // neighbor mz
                thrust::make_permutation_iterator(_spectrum.fI->begin(), seqindex.begin()), 
                thrust::make_permutation_iterator(_miUsed.begin(), neighbors.begin() ),
                dscores.begin())), 
        thrust::make_zip_iterator(
            thrust::make_tuple(
                neighbors.end(),
                thrust::make_permutation_iterator(_spectrum.iM->begin(), seqindex.end()), 
                thrust::make_permutation_iterator(_spectrum.fI->begin(), seqindex.end()), 
                thrust::make_permutation_iterator(_miUsed.begin(), neighbors.end()),
                dscores.end())), 
            dot_neighbor_functor());
    // look right 
    // 1+lCount thing reproduces k-score behavior when actual ion does not match but neighbors do
    thrust::transform(lcounts.begin(),lcounts.end(),lcounts.begin(),incrementby(1));
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(seqindex.begin(), lcounts.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(seqindex.end(), lcounts.end())), add_B_to_A_functor());
    thrust::transform(neighbors.begin(),neighbors.end(),neighbors.begin(),incrementby(2));
    thrust::for_each(
        thrust::make_zip_iterator( 
            thrust::make_tuple(
                neighbors.begin(), 
                thrust::make_permutation_iterator(_spectrum.iM->begin(), seqindex.begin()), // neighbor mz
                thrust::make_permutation_iterator(_spectrum.fI->begin(), seqindex.begin()), 
                thrust::make_permutation_iterator(_miUsed.begin(), neighbors.begin()),
                dscores.begin())), // neighbor contributor flag
        thrust::make_zip_iterator(
            thrust::make_tuple(
                neighbors.end(),
                thrust::make_permutation_iterator(_spectrum.iM->begin(), seqindex.end()), 
                thrust::make_permutation_iterator(_spectrum.fI->begin(), seqindex.end()), 
                thrust::make_permutation_iterator(_miUsed.begin(), neighbors.end()),
                dscores.end())), 
            dot_neighbor_functor());
    // and get total score
    float dScore = thrust::reduce(dscores.begin(),dscores.end());
    return (dScore);
}

