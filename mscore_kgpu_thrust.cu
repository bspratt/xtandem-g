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
#include <thrust/adjacent_difference.h>
#include <map>

#ifdef HAVE_MULTINODE_TANDEM // support for Hadoop and/or MPI?
#undef HAVE_MULTINODE_TANDEM // need to omit boost etc for NVCC's benefit
#endif
using namespace std; // this is evil but tandem codebase assumes
typedef map<unsigned long,double> SMap;
#include "mscore_kgpu_thrust.h"

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

#define ZEROVEC(f) thrust::fill((f).begin(),(f).end(),0)  // TODO - cudaMemSet?

void mscore_kgpu_thrust_fvec_clear(thrust::device_vector<float> *vec) {
    ZEROVEC(*vec);
}
void mscore_kgpu_thrust_fvec_kill(thrust::device_vector<float> *vec) {
    delete vec;
}

void mscore_kgpu_thrust_ivec_kill(thrust::device_vector<int> *vec) {
    delete vec;
}

#define MY_CUDA_STREAM 0

thrust::device_vector<int> *mscore_kgpu_thrust_host_to_device_copy(const pinned_host_vector_int_t &host_src,thrust::device_vector<int> *device_dest) {
    if (!device_dest) {
        device_dest = new thrust::device_vector<int>(host_src.size());
    } else {
        device_dest->resize(host_src.size());
    }
    cudaMemcpyAsync(thrust::raw_pointer_cast(device_dest->data()),host_src.data(),host_src.size()*sizeof(int),cudaMemcpyHostToDevice,MY_CUDA_STREAM);
    return device_dest;
}

void mscore_kgpu_thrust_host_to_device_copy_float(const pinned_host_vector_float_t &host_src,thrust::device_vector<float> &device_dest) {
    device_dest.resize(host_src.size());
    cudaMemcpyAsync(thrust::raw_pointer_cast(device_dest.data()),host_src.data(),host_src.size()*sizeof(float),cudaMemcpyHostToDevice,MY_CUDA_STREAM);
}

void mscore_kgpu_thrust_device_to_host_copy_float(const thrust::device_vector<float> &device_src,pinned_host_vector_float_t &host_dest) {
    host_dest.resize(device_src.size());
    cudaMemcpyAsync(thrust::raw_pointer_cast(host_dest.data()),thrust::raw_pointer_cast(device_src.data()),host_dest.size()*sizeof(float),cudaMemcpyDeviceToHost,MY_CUDA_STREAM);
}

void mscore_kgpu_thrust_device_to_host_copy_int(const thrust::device_vector<int> &device_src,pinned_host_vector_int_t &host_dest) {
    host_dest.resize(device_src.size());
    cudaMemcpyAsync(thrust::raw_pointer_cast(host_dest.data()),thrust::raw_pointer_cast(device_src.data()),host_dest.size()*sizeof(int),cudaMemcpyDeviceToHost,MY_CUDA_STREAM);
}

#define HOSTCODE // __host__

// define transformation f(x) -> x^2
struct square
{
    HOSTCODE __device__
        float operator()(float x)
    {
        return x * x;
    }
};

struct divideby
{
    const float m_factor;
    HOSTCODE __device__
        divideby(float divisor) :  m_factor (1.0f/divisor) {}
    HOSTCODE __device__
        float operator()(float m)
    {
        return (m*m_factor);
    }
};

struct multiplyby
{
    const float m_factor;
    HOSTCODE __device__
        multiplyby(float factor) : m_factor(factor) {}
    HOSTCODE __device__
        float operator()(float m)
    {
        return (m*m_factor);
    }
};

struct xform_multiplyby  : public thrust::unary_function<int,int>
{
    const int m_factor;
    HOSTCODE __device__
        xform_multiplyby(int factor) : m_factor(factor) {}
    HOSTCODE __device__
        int operator()(int m)
    {
        return (m*m_factor);
    }
};


struct xform_incrementby  : public thrust::unary_function<int,int>
{
    const int m_incr;
    HOSTCODE __device__
        xform_incrementby(int incr) :  m_incr (incr) {}
    HOSTCODE __device__
        int operator()(int m)
    {
        return (m+m_incr);
    }
};

struct xform_modby  : public thrust::unary_function<int,int>
{
    const int m_mod;
    HOSTCODE __device__
        xform_modby(int m) :  m_mod (m) {}
    HOSTCODE __device__
        int operator()(int m)
    {
        return (m%m_mod);
    }
};



/*
* m/z binning
*/
struct imass_functor
{
    const float m_IsotopeCorrection;
    const int m_binOffset;
    HOSTCODE __device__
        imass_functor(float IsotopeCorrection, int binOffset) : 
    m_IsotopeCorrection((float)(1.0/IsotopeCorrection)),m_binOffset(binOffset) {}
    template <typename Tuple>
    HOSTCODE __device__
        void operator()(Tuple m) // m[0]=raw mz m[1]=binned mz
    {
        thrust::get<1>(m) =  ((int)((thrust::get<0>(m)*m_IsotopeCorrection) + 0.5f))-m_binOffset;
    }
};


/*
* mixrange
* fI[i] = fi_scattered[i] - (sums[i+101]-sums[i])/101
* <0> = <1> - (<2>-<3>)/101
*/
struct mixrange_functor
{
    template <typename Tuple>
    HOSTCODE __device__
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
    HOSTCODE __device__
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

struct a_equals_b_minus_c
{
    template <typename Tuple>
    HOSTCODE __device__
        void operator()(Tuple m) 
    {
        thrust::get<0>(m) = thrust::get<1>(m)-thrust::get<2>(m);
    }
};

struct nonpositive_functor
{
    template <typename Tuple>
    HOSTCODE __device__
        bool operator()(Tuple m) 
    {
        return (thrust::get<1>(m) <= 0);
    }
};

// define transformation f(x) -> sqrt(x)
struct squareroot
{
    HOSTCODE __device__
        float operator()(float x)
    {
        return sqrt(x);
    }
};

CUDA_TIMER_DECL(tScore);

void mscore_kgpu_thrust_score(
    // m/z-intensity pairs (input)
    const thrust::device_vector<float>::iterator itfM,
    const thrust::device_vector<float>::iterator itfI,    
    size_t length, // iterator end distance
    double dIsotopeCorrection, // for m/z binning (input)
    int iWindowCount, // (input)
    int endMassMax, // max acceptable binned m/z (input)
    int &m_maxEnd,  // max encountered binned m/z (output)
    vmiTypeGPU &binned)   // binned intensities and m/z's (output)
{
CUDA_TIMER_START(tScore)
    if (!length) {
        return;
    }
    thrust::device_vector<int> &iM = *binned.iM;
    thrust::device_vector<float> &fI = *binned.fI;
    thrust::copy_n(itfI,length,fI.begin());

    // figure the shift needed for lowest binned mass and half smoothing window
    int startMass = (int)((*itfM*dIsotopeCorrection) + 0.5f);
    const int binOffset = startMass-50;
    endMassMax -= binOffset; // shift to local coords
    startMass -= binOffset;

    // bin the m/z's
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(itfM, iM.begin())), 
                     thrust::make_zip_iterator(thrust::make_tuple(itfM+length, iM.end())), 
                     imass_functor(dIsotopeCorrection,binOffset));

    // Screen peeks on upper end.  TODO faster way? thrust::lowerbound is slower than this
    iveciter itM = iM.begin();
    iveciter itEnd = iM.end();
    int endMass = itEnd[-1];
    int trimEnd=0; 
    while (itM != itEnd && endMass >= endMassMax) {
        itEnd--;
        trimEnd++;
        endMass = itEnd[-1];
    }

    if (itM == itEnd)  // No peaks left.
    {
        return;
    }

    // if any bin has more than one occupant, choose the one with higher intensity
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(itM, itM+1, fI.begin(), fI.begin()+1)), 
                     thrust::make_zip_iterator(thrust::make_tuple(itEnd-1, itEnd, fI.end()-(trimEnd+1), fI.end()-trimEnd)), 
                     imass_collision_functor());

    // take sqrt of intensities
    thrust::transform(fI.begin(),fI.end()-trimEnd,fI.begin(),squareroot());

    fveciter itI = fI.begin();

    // for every pair fI,iM set templookup[iM]=fI
    thrust::device_vector<float> fI_scattered(endMassMax+50); // +50 for smoothing window
    thrust::scatter(fI.begin(),fI.end()-trimEnd,iM.begin(),fI_scattered.begin());
    if ((endMass+binOffset)+50 > m_maxEnd)
        m_maxEnd = (endMass+binOffset)+50;

    float fMaxI = *thrust::max_element(fI.begin(),fI.end()-trimEnd);
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
    // TODO make this loop into a single call
    /*
     * Process input spectrum for dot product - split windows
     */
    for (int i = 0; i < iWindowCount; i++) {

        int iStart = startMass + i*iWindowSize;

        /*
         * Get maximum intensity within window
         */
        float fMaxWindowI = *thrust::max_element(fI_scattered.begin()+iStart,fI_scattered.begin()+(iStart+iWindowSize));

        if (fMaxWindowI > 0.0 && fMaxWindowI > fMinCutoff) {
            /*
             * Normalize within window
             */
            thrust::transform(fI_scattered.begin()+iStart,fI_scattered.begin()+(iStart + iWindowSize),
                fI_scattered.begin()+iStart,multiplyby(fMaxI / fMaxWindowI));
        }
    }

 /*
    * Reduce intensity and make unit vector by dividing
    * every point by sqrt(sum(x^2))
    */
    double dSpectrumArea = sqrt(thrust::transform_reduce(fI_scattered.begin()+startMass, fI_scattered.begin()+(endMass+1), square(), 0.0f, thrust::plus<float>()));
    thrust::transform(fI_scattered.begin()+startMass,fI_scattered.begin()+(endMass+1),fI_scattered.begin()+startMass,divideby(dSpectrumArea));

    /*
     * Perform mix-range modification to input spectrum
     */
    thrust::device_vector<float> sums(fI_scattered.size()+101);
    thrust::copy(fI_scattered.begin(),fI_scattered.end(),sums.begin()+50);
    thrust::exclusive_scan(sums.begin()+50, sums.end(), sums.begin()+50);
    // now sums[i+101]-sums[i] = sum(fI_scattered[i-50:i+50])
    fI.resize(fI_scattered.size());
    // fI[i] = fI_scattered[i] - (sums[i+101]-sums[i])/101
    thrust::for_each(thrust::make_zip_iterator(
                     thrust::make_tuple( fI.begin(), fI_scattered.begin(), sums.begin()+101, sums.begin())),
                     thrust::make_zip_iterator(
                     thrust::make_tuple( fI.end(), fI_scattered.end(), sums.end(),sums.end()-101)), 
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
    binned.iM_max = iM.back();
CUDA_TIMER_STOP(tScore)
}


/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses).
 */

struct dot_functor
{
    HOSTCODE __device__
        dot_functor() {}
    template <typename Tuple>
    HOSTCODE __device__
        void operator()(Tuple m) // 0=intensity, 1=sequence (0, .5, or 1) 2=used, 3=lcount
    {
        if (thrust::get<1>(m)) { // ion match?
            float contrib = thrust::get<1>(m)*thrust::get<0>(m); // intensity * ion contrib (1.0 for match, .5 for neighbor)
            if (thrust::get<2>(m) < contrib) {
                thrust::get<3>(m) += (1.0==thrust::get<1>(m)); // lcount++ if direct ion match
                thrust::get<2>(m) = contrib; // m_used
            }
        }
    }
};

#define MAX_SEQLEN 200
thrust::device_vector<float> ones;
thrust::device_vector<float> halves;
thrust::device_vector<float> seq_hits;
int seq_hits_len=0; // seq_hits is a concatenation of lists each seq_hits_len long
thrust::device_vector<float> dScoresDev;
thrust::device_vector<int> lCountsDev;
thrust::device_vector<int> lcounts;
thrust::device_vector<int> lcountsscan;
thrust::device_vector<float> miUsed;
thrust::device_vector<float> scatteredCopies;

void mscore_kgpu_thrust_init() {
  ones.resize(MAX_SEQLEN);
  thrust::fill(ones.begin(), ones.end(), 1.0);
  halves.resize(MAX_SEQLEN);
  thrust::fill(halves.begin(), halves.end(), 0.5);
}

// perform dot on current spectrum and all sequence variations at one go
void mscore_kgpu_thrust_dot(pinned_host_vector_int_t &lCountsResult,
    pinned_host_vector_float_t &dScoresResult,
    const std::vector<const vmiTypeGPU *> &spectra,
    const std::vector<int> &sequenceCountAfterEachScorePreload,
    const thrust::device_vector<int> &cached_sequences, 
    const std::vector<int> &vsequence_index)
{
//puts("mscore_kgpu_thrust_dot");
    cudaStreamSynchronize(MY_CUDA_STREAM); // wait for any preperatory memcpy to complete
    int lastEnd = 0;
    size_t nSeqTotal = vsequence_index.size()-1;
    dScoresDev.resize(nSeqTotal);
    lCountsDev.resize(nSeqTotal);
    for (int spec=0;spec<spectra.size();lastEnd=sequenceCountAfterEachScorePreload[spec++]) {
        const int *sequence_index = &vsequence_index[lastEnd];
        int nSeq = (sequenceCountAfterEachScorePreload[spec]-lastEnd);
        const vmiTypeGPU &spectrum = *spectra[spec];
        int seq_hits_len = spectrum.iM_max+1; // sequence is trimmed to max spectrum mz
        if (lcounts.size() < seq_hits_len*nSeq) {
            lcounts.resize((seq_hits_len*nSeq)+1); // need a -1th element for sum diff later
            lcountsscan.resize((seq_hits_len*nSeq)+1); // need a -1th element for sum diff later
            seq_hits.resize(seq_hits_len*nSeq);
            miUsed.resize((seq_hits_len*nSeq)+1); // need a -1th element for sum diff later
            scatteredCopies.resize(seq_hits_len*nSeq);
        }
        ZEROVEC(seq_hits);
        // set up the .5 and 1.0 sequence hit multipliers
        for (int s=0;s<nSeq;s++ ) {
            int dist = sequence_index[s+1]-sequence_index[s];
            int seq_hits_offset = s*seq_hits_len;
            thrust::device_vector<int>::const_iterator seq_begin = cached_sequences.begin()+sequence_index[s];
            // note the .5 contributions for ions to left and right of actual ions
            thrust::scatter(halves.begin(),halves.begin()+dist,thrust::make_transform_iterator(seq_begin,xform_incrementby(-1)),seq_hits.begin()+seq_hits_offset);
            thrust::scatter(halves.begin(),halves.begin()+dist,thrust::make_transform_iterator(seq_begin,xform_incrementby(1)),seq_hits.begin()+seq_hits_offset);
            for (int ss=s+1;ss<nSeq;ss++ ) { // and make it an underlayment for following sequences
                int seq_hits_offset = ss*seq_hits_len;
                thrust::scatter(halves.begin(),halves.begin()+dist,thrust::make_transform_iterator(seq_begin,xform_incrementby(-1)),seq_hits.begin()+seq_hits_offset);
                thrust::scatter(halves.begin(),halves.begin()+dist,thrust::make_transform_iterator(seq_begin,xform_incrementby(1)),seq_hits.begin()+seq_hits_offset);
            }
        }
        for (int s=0;s<nSeq;s++ ) {
            int dist = sequence_index[s+1]-sequence_index[s];
            int seq_hits_offset = s*seq_hits_len;
            thrust::device_vector<int>::const_iterator seq_begin = cached_sequences.begin()+sequence_index[s];
            // note the 1.0 contribution of actual ions
            thrust::scatter(ones.begin(),ones.begin()+dist,seq_begin,seq_hits.begin()+seq_hits_offset);
            for (int ss=s+1;ss<nSeq;ss++ ) { // and make it an underlayment for following sequences
                int seq_hits_offset = ss*seq_hits_len;
                thrust::scatter(ones.begin(),ones.begin()+dist,seq_begin,seq_hits.begin()+seq_hits_offset);
            }
        }
        // now lay out a string of spectrum copies so we can do all sequences in one shot
        ZEROVEC(scatteredCopies);
        for (int s=0;s<nSeq;s++ ) {
            thrust::scatter(spectrum.fI->begin(),spectrum.fI->end(),spectrum.iM->begin(),scatteredCopies.begin()+(s*seq_hits_len));
        }


        ZEROVEC(lcounts);
        ZEROVEC(miUsed);
        // now find the hits
        thrust::for_each(
            thrust::make_zip_iterator( 
            thrust::make_tuple(
            scatteredCopies.begin(), 
            seq_hits.begin(), 
            miUsed.begin()+1, 
            lcounts.begin()+1)), 
            thrust::make_zip_iterator(
            thrust::make_tuple(
            scatteredCopies.end(), 
            seq_hits.end(), 
            miUsed.end(),
            lcounts.end())), 
            dot_functor());
        // and get total score 
        thrust::inclusive_scan(lcounts.begin(),lcounts.end(),lcountsscan.begin()); // TODO maybe count_if?
        thrust::for_each(
            thrust::make_zip_iterator( thrust::make_tuple(
            lcounts.begin(),
            thrust::make_permutation_iterator(lcountsscan.begin(),
            thrust::make_transform_iterator(thrust::make_counting_iterator(0),xform_multiplyby(seq_hits_len))+1),
            thrust::make_permutation_iterator(lcountsscan.begin(),
            thrust::make_transform_iterator(thrust::make_counting_iterator(0),xform_multiplyby(seq_hits_len))))),
            thrust::make_zip_iterator( thrust::make_tuple(
            lcounts.begin()+nSeq,
            thrust::make_permutation_iterator(lcountsscan.begin(),
            thrust::make_transform_iterator(thrust::make_counting_iterator(0),xform_multiplyby(seq_hits_len))+nSeq+1),
            thrust::make_permutation_iterator(lcountsscan.begin(),
            thrust::make_transform_iterator(thrust::make_counting_iterator(0),xform_multiplyby(seq_hits_len))+nSeq))),
            a_equals_b_minus_c());
        // lcounts now contains count of 1.0 hits for each seq, but we want just the increment
        thrust::adjacent_difference(lcounts.begin(),lcounts.begin()+nSeq,lCountsDev.begin()+lastEnd);

        // now find sum of miUsed for each seq
        thrust::inclusive_scan(miUsed.begin()+1,miUsed.end(),miUsed.begin()+1);
        thrust::for_each(
            thrust::make_zip_iterator( thrust::make_tuple(
            miUsed.begin()+1,
            thrust::make_permutation_iterator(miUsed.begin(),
            thrust::make_transform_iterator(thrust::make_counting_iterator(0),xform_multiplyby(seq_hits_len))+1),
            thrust::make_permutation_iterator(miUsed.begin(),
            thrust::make_transform_iterator(thrust::make_counting_iterator(0),xform_multiplyby(seq_hits_len))))),
            thrust::make_zip_iterator( thrust::make_tuple(
            miUsed.begin()+1+nSeq,
            thrust::make_permutation_iterator(miUsed.begin(),
            thrust::make_transform_iterator(thrust::make_counting_iterator(0),xform_multiplyby(seq_hits_len))+nSeq+1),
            thrust::make_permutation_iterator(miUsed.begin(),
            thrust::make_transform_iterator(thrust::make_counting_iterator(0),xform_multiplyby(seq_hits_len))+nSeq))),
            a_equals_b_minus_c());
        // miUsed[1:n] now contains sum of hits for each seq, but we want just the increment
        thrust::adjacent_difference(miUsed.begin()+1,miUsed.begin()+1+nSeq,dScoresDev.begin()+lastEnd);
    }
    mscore_kgpu_thrust_device_to_host_copy_int(lCountsDev,lCountsResult); // copy back to host memory
    mscore_kgpu_thrust_device_to_host_copy_float(dScoresDev,dScoresResult); // copy back to host memory
    cudaStreamSynchronize(MY_CUDA_STREAM); // wait for memcpy to complete
    return;
}

