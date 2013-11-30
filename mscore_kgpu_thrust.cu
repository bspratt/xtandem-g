// -*- mode: c++ -*-

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
#define CLEARVEC(f) (f).resize(0);(f).shrink_to_fit();
static void clear_largebuffers(); // defined below

template<typename vectT> bool RESIZEVEC(vectT &f,size_t sz,bool recursionOK=true) {
    try {
        f.resize(sz);
    }
    catch(exception &e) {
        // fail - try to hose out some larger buffers 
        if (recursionOK && !largeBufferLockEngaged)
        {
            clear_largebuffers();
            RESIZEVEC(f,sz,false);
        }
        return false;
    }
    return true;
}

void mscore_kgpu_thrust_fvec_clear(thrust::device_vector<float> *vec) {
    ZEROVEC(*vec);
}
void mscore_kgpu_thrust_fvec_kill(thrust::device_vector<float> *vec) {
    delete vec;
}

void mscore_kgpu_thrust_ivec_kill(thrust::device_vector<int> *vec) {
    delete vec;
}

// convert a linear index to a row index
template <typename T>
struct linear_index_to_row_index : public thrust::unary_function<T,T>
{
    T C; // number of columns
    
    __host__ __device__
    linear_index_to_row_index(T _C) : C(_C) {}

    __host__ __device__
    T operator()(T i)
    {
        return i / C;
    }
};

#define MY_CUDA_STREAM 0

thrust::device_vector<int> *mscore_kgpu_thrust_host_to_device_copy(
			      const pinned_host_vector_int_t &host_src,
			      thrust::device_vector<int> *device_dest) {
    if (!device_dest) {
        device_dest = new thrust::device_vector<int>();
    }
    RESIZEVEC(*device_dest,host_src.size());
    cudaMemcpyAsync(thrust::raw_pointer_cast(device_dest->data()),
		    host_src.data(),
		    host_src.size()*sizeof(int),
		    cudaMemcpyHostToDevice,
		    MY_CUDA_STREAM);
    return device_dest;
}

void mscore_kgpu_thrust_host_to_device_copy_float(
       const pinned_host_vector_float_t &host_src,
       thrust::device_vector<float> &device_dest) {
    RESIZEVEC(device_dest,host_src.size());
    cudaMemcpyAsync(thrust::raw_pointer_cast(device_dest.data()),
		    host_src.data(),
		    host_src.size()*sizeof(float),
		    cudaMemcpyHostToDevice,
		    MY_CUDA_STREAM);
}

void mscore_kgpu_thrust_device_to_host_copy_float(
	const thrust::device_vector<float> &device_src,
	pinned_host_vector_float_t &host_dest) {
    RESIZEVEC(host_dest,device_src.size());
    cudaMemcpyAsync(thrust::raw_pointer_cast(host_dest.data()),
		    thrust::raw_pointer_cast(device_src.data()),
		    device_src.size()*sizeof(float),
		    cudaMemcpyDeviceToHost,
		    MY_CUDA_STREAM);
}

void mscore_kgpu_thrust_device_to_host_copy_int(
	const thrust::device_vector<int> &device_src,
	pinned_host_vector_int_t &host_dest) {
    RESIZEVEC(host_dest,device_src.size());
    cudaMemcpyAsync(thrust::raw_pointer_cast(host_dest.data()),
		    thrust::raw_pointer_cast(device_src.data()),
		    device_src.size()*sizeof(int),
		    cudaMemcpyDeviceToHost,
		    MY_CUDA_STREAM);
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
        thrust::get<1>(m) =
	  ((int)((thrust::get<0>(m)*m_IsotopeCorrection) + 0.5f)) - m_binOffset;
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
        thrust::get<0>(m) =
	  thrust::get<1>(m) -  ((thrust::get<2>(m)-thrust::get<3>(m))/101.0f);
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
    thrust::for_each(thrust::make_zip_iterator(
		       thrust::make_tuple(itfM, iM.begin())),
		     thrust::make_zip_iterator(
			thrust::make_tuple(itfM+length, iM.end())),
                     imass_functor(dIsotopeCorrection,binOffset));

    // Screen peeks on upper end.  
    // TODO faster way? thrust::lowerbound is slower than this
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

    // if any bin has more than one occupant, choose the one with
    // higher intensity
    thrust::for_each(
	thrust::make_zip_iterator(thrust::make_tuple(itM, 
						     itM+1, 
						     fI.begin(), 
						     fI.begin()+1)),
	thrust::make_zip_iterator(thrust::make_tuple(itEnd-1,
						     itEnd,
						     fI.end()-(trimEnd+1),
						     fI.end()-trimEnd)), 
	imass_collision_functor());

    // take sqrt of intensities
    thrust::transform(fI.begin(),fI.end()-trimEnd,fI.begin(),squareroot());

    fveciter itI = fI.begin();

    // for every pair fI,iM set templookup[iM]=fI
    thrust::device_vector<float> fI_scattered(endMassMax+50); // +50
							      // for
							      // smoothing
							      // window
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
        float fMaxWindowI = *thrust::max_element(
	  fI_scattered.begin()+iStart,
	  fI_scattered.begin()+(iStart+iWindowSize));

        if (fMaxWindowI > 0.0 && fMaxWindowI > fMinCutoff) {
            /*
             * Normalize within window
             */
            thrust::transform(
	      fI_scattered.begin()+iStart,
	      fI_scattered.begin()+(iStart + iWindowSize),
	      fI_scattered.begin()+iStart,multiplyby(fMaxI / fMaxWindowI));
        }
    }

 /*
    * Reduce intensity and make unit vector by dividing
    * every point by sqrt(sum(x^2))
    */
    double dSpectrumArea = sqrt(thrust::transform_reduce(
				  fI_scattered.begin()+startMass, 
				  fI_scattered.begin()+(endMass+1), 
				  square(), 
				  0.0f, 
				  thrust::plus<float>()));
    thrust::transform(fI_scattered.begin()+startMass,
		      fI_scattered.begin()+(endMass+1),
		      fI_scattered.begin()+startMass,
		      divideby(dSpectrumArea));

    /*
     * Perform mix-range modification to input spectrum
     */
    thrust::device_vector<float> sums(fI_scattered.size()+101);
    thrust::copy(fI_scattered.begin(),fI_scattered.end(),sums.begin()+50);
    thrust::exclusive_scan(sums.begin()+50, sums.end(), sums.begin()+50);
    // now sums[i+101]-sums[i] = sum(fI_scattered[i-50:i+50])
    RESIZEVEC(fI,fI_scattered.size());
    // fI[i] = fI_scattered[i] - (sums[i+101]-sums[i])/101
    thrust::for_each(thrust::make_zip_iterator(
                     thrust::make_tuple(fI.begin(), 
					fI_scattered.begin(), 
					sums.begin()+101, 
					sums.begin())),
                     thrust::make_zip_iterator(
		       thrust::make_tuple(
			 fI.end(), 
			 fI_scattered.end(), 
			 sums.end(),
			 sums.end()-101)), 
                     mixrange_functor());

    // now remove any non-positive results
    RESIZEVEC(iM,fI.size());
    thrust::sequence(iM.begin(),iM.end(),binOffset,1); // shift back
						       // to host
						       // coords
     thrust::zip_iterator<thrust::tuple<iveciter,fveciter> > new_end = 
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
        void operator()(Tuple m) // 0=intensity, 1=sequence (0, .5, or
				 // 1) 2=used, 3=lcount
    {
        if (thrust::get<1>(m)) { // ion match?

	    // intensity * ion contrib (1.0 for match, .5 for neighbor)
            float contrib = thrust::get<1>(m)*thrust::get<0>(m); 
            
	    if (thrust::get<2>(m) < contrib) {

		// lcount++ if direct ion match
                thrust::get<3>(m) += (1.0==thrust::get<1>(m));

                thrust::get<2>(m) = contrib; // m_used
            }
        }
    }
};

#define MAX_SEQLEN 200
thrust::device_vector<float> ones;
thrust::device_vector<float> halves;
thrust::device_vector<float> seq_hits;
// seq_hits is a concatenation of lists each seq_hits_len long
thrust::device_vector<float> dScoresDev;
thrust::device_vector<int> lcounts;
thrust::device_vector<int> lCountsDev;
thrust::device_vector<int> lcountsTmp;
thrust::device_vector<float> miUsed;
thrust::device_vector<float> miUsedTmp;
thrust::device_vector<float> scatteredCopies;
thrust::device_vector<int> row_indices;

//
// helps with the problem of large allocs
// preventing smaller ones - we can break
// those big ones up if needed
//
static bool largeBufferLockEngaged = false;
class LargeBufferLock {
    LargeBufferLock() {
        largeBufferLockEngaged = true;
    }
    ~LargeBufferLock() {
        largeBufferLockEngaged = false;
    }
};

static void clear_largebuffers() { 
    // clear out our larger buffers if possible
    if (!largeBufferLockEngaged) {
        CLEARVEC(lcountsTmp);
        CLEARVEC(seq_hits);
        CLEARVEC(miUsedTmp);
        CLEARVEC(scatteredCopies);
    }
}

void mscore_kgpu_thrust_init(size_t available_memsize) {
  RESIZEVEC(ones,MAX_SEQLEN);
  thrust::fill(ones.begin(), ones.end(), 1.0);
  RESIZEVEC(halves,MAX_SEQLEN);
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
    //cudaStreamSynchronize(MY_CUDA_STREAM); // wait for any preperatory
					   // memcpy to complete
    int lastEnd = 0;
    size_t nSeqTotal = vsequence_index.size()-1;
    RESIZEVEC(dScoresDev,nSeqTotal);
    RESIZEVEC(lCountsDev,nSeqTotal);
    RESIZEVEC(lcounts,nSeqTotal); 
    RESIZEVEC(miUsed,nSeqTotal);
    RESIZEVEC(row_indices,nSeqTotal);
    int seq_hits_len = 0;
    for (int sp=(int)spectra.size();sp--;) {
        seq_hits_len = max(seq_hits_len,spectra[sp]->iM_max+1); 
    }
    seq_hits_len = 32*(1+(seq_hits_len/32)); // nice even boundaries
    int len = seq_hits_len*nSeqTotal;
    int nSpectraPerChunk = (int)spectra.size();
    LargeBufferLock lock();
    if ((int)lcountsTmp.size() < len) {
//std::cout<<len<<"\n";fflush(stdout);
      while (!(RESIZEVEC(lcountsTmp,len) && 
           RESIZEVEC(seq_hits,len) &&
           RESIZEVEC(miUsedTmp,len) && 
           RESIZEVEC(scatteredCopies,len)))
      {
          // too big for memory, break task into smaller chunks
          nSpectraPerChunk /= 2;
          nSpectraPerChunk += (int)(spectra.size()%2); // 999->500+499 instead of 499+499+1
          int maxSeqPerChunk=0;
          for (int spec=0;spec<(int)spectra.size();) {
              // determine worst case memory use
              int firstSeqInChunk = spec?sequenceCountAfterEachScorePreload[spec-1]:0;
              int nextspec = min((int)spectra.size(),spec+nSpectraPerChunk);
              int lastSeqInChunk = sequenceCountAfterEachScorePreload[nextspec-1];
              maxSeqPerChunk = max(maxSeqPerChunk,lastSeqInChunk-firstSeqInChunk);
              spec = nextspec;
          }
          len = seq_hits_len*maxSeqPerChunk;
      }
    }

    // in case of looping, use this to incrementally copy back into overall result space
    int devOffset = 0; 

    for (int firstSpectraInChunk=0;firstSpectraInChunk<(int)spectra.size();firstSpectraInChunk+=nSpectraPerChunk)
    {
        // normally this loop just hits once, but written this way to subdivide
        // larger than memory items

        int firstSpectraNotInChunk = min((int)spectra.size(),firstSpectraInChunk+nSpectraPerChunk);
        ZEROVEC(seq_hits);
        ZEROVEC(scatteredCopies);
        ZEROVEC(lcountsTmp);
        ZEROVEC(miUsedTmp);

        int row = 0;
        int rowB = 0;
        for (int spec=firstSpectraInChunk;
	     spec<firstSpectraNotInChunk;
	     lastEnd=sequenceCountAfterEachScorePreload[spec++]) {
             int rowA = row;
            const int *sequence_index = &vsequence_index[lastEnd];
            int nSeq = (sequenceCountAfterEachScorePreload[spec]-lastEnd);
            const vmiTypeGPU &spectrum = *spectra[spec];
            // set up the .5 and 1.0 sequence hit multipliers
            for (int s=0;s<nSeq;s++ ) {
                int dist = sequence_index[s+1]-sequence_index[s];
                int seq_hits_offset = row++*seq_hits_len;
                thrust::device_vector<int>::const_iterator seq_begin = 
	          cached_sequences.begin()+sequence_index[s];
            
	        // note the .5 contributions for ions to left and right of actual ions
                thrust::scatter(
	          halves.begin(),
	          halves.begin()+dist,
	          thrust::make_transform_iterator(seq_begin,
					          xform_incrementby(-1)),
	          seq_hits.begin()+seq_hits_offset);
	    
                thrust::scatter(
	          halves.begin(),
	          halves.begin()+dist,
	          thrust::make_transform_iterator(seq_begin,
					          xform_incrementby(1)),
	          seq_hits.begin()+seq_hits_offset);

                for (int ss=s+1;ss<nSeq;ss++ ) { // and make it an
					         // underlayment for
					         // following sequences
                    int seq_hits_offset = (rowA+ss)*seq_hits_len;
                    thrust::scatter(halves.begin(),
				    halves.begin()+dist,
				    thrust::make_transform_iterator(
				      seq_begin, xform_incrementby(-1)),
				    seq_hits.begin()+seq_hits_offset);

                    thrust::scatter(halves.begin(),
				    halves.begin()+dist,
				    thrust::make_transform_iterator(
				      seq_begin,xform_incrementby(1)),
				    seq_hits.begin()+seq_hits_offset);
                }
            }
            for (int s=0;s<nSeq;s++ ) {
                int dist = sequence_index[s+1]-sequence_index[s];
                int seq_hits_offset = (rowA+s)*seq_hits_len;
                thrust::device_vector<int>::const_iterator seq_begin = 
	          cached_sequences.begin()+sequence_index[s];

                // note the 1.0 contribution of actual ions
                thrust::scatter(ones.begin(),
			        ones.begin()+dist,
			        seq_begin,
			        seq_hits.begin()+seq_hits_offset);
                for (int ss=s+1;ss<nSeq;ss++ ) { // and make it an
					         // underlayment for
					         // following sequences
                    int seq_hits_offset = (rowA+ss)*seq_hits_len;
                    thrust::scatter(ones.begin(),
				    ones.begin()+dist,
				    seq_begin,
				    seq_hits.begin()+seq_hits_offset);
                }
            }
            // now lay out a string of spectrum copies so we can do all sequences in one shot
            for (int s=0;s<nSeq;s++ ) {
                thrust::scatter(spectrum.fI->begin(),
			        spectrum.fI->end(),
			        spectrum.iM->begin(),
			        scatteredCopies.begin()+(rowB++*seq_hits_len));
            }
        } // end for each spectra
        // now find the hits
        thrust::for_each(
            thrust::make_zip_iterator( 
            thrust::make_tuple(
            scatteredCopies.begin(), 
            seq_hits.begin(), 
            miUsedTmp.begin()+1, 
            lcountsTmp.begin()+1)), 
            thrust::make_zip_iterator(
            thrust::make_tuple(
            scatteredCopies.end(), 
            seq_hits.end(), 
            miUsedTmp.end(),
            lcountsTmp.end())), 
            dot_functor());
#ifdef ___DEBUG
STDVECT(int,lcountsTmpLocal,lcountsTmp);
printf("lcountsTmpLocal ");for (int nnn=0;nnn<lcountsTmpLocal.size();nnn++)if(lcountsTmpLocal[nnn])printf("%d,%d ",nnn,lcountsTmpLocal[nnn]);printf("\n\n");fflush(stdout);
#endif
        // and get total score 
        // compute row sums by summing values with equal row indices
        int firstSeqInChunk = firstSpectraInChunk?sequenceCountAfterEachScorePreload[firstSpectraInChunk-1]:0;
        int nSeqInChunk = sequenceCountAfterEachScorePreload[firstSpectraNotInChunk-1] - firstSeqInChunk;
        thrust::reduce_by_key(thrust::make_transform_iterator(thrust::counting_iterator<int>(0), linear_index_to_row_index<int>(seq_hits_len)),
                              thrust::make_transform_iterator(thrust::counting_iterator<int>(0), linear_index_to_row_index<int>(seq_hits_len)) + (seq_hits_len*(int)nSeqInChunk),
                              lcountsTmp.begin(),
                              row_indices.begin(),
                              lcounts.begin()+firstSeqInChunk,
                              thrust::equal_to<int>(),
                              thrust::plus<int>());
#ifdef ___DEBUG
STDVECT(int,lcountss,lcounts);
printf("lcounts ");for (int nnn=0;nnn<lcountss.size();nnn++)printf("%d,%d ",nnn,lcountss[nnn]);printf("\n\n");fflush(stdout);
#endif

        // now find sum of miUsedTmp for each seq
        // compute row sums by summing values with equal row indices
        thrust::reduce_by_key(thrust::make_transform_iterator(thrust::counting_iterator<int>(0), linear_index_to_row_index<int>(seq_hits_len)),
                             thrust::make_transform_iterator(thrust::counting_iterator<int>(0), linear_index_to_row_index<int>(seq_hits_len)) + (seq_hits_len*(int)nSeqInChunk),
                             miUsedTmp.begin()+1,
                             row_indices.begin(),
                             miUsed.begin()+firstSeqInChunk,
                             thrust::equal_to<int>(),
                             thrust::plus<float>());
#ifdef ___DEBUG
STDVECT(float,miUsedTmpLocal,miUsed);
printf("miUsed ");for (int nnn=0;nnn<miUsedTmpLocal.size();nnn++)if(miUsedTmpLocal[nnn])printf("%d,%f ",nnn,miUsedTmpLocal[nnn]);printf("\n\n");fflush(stdout);
#endif

        // miUsed[0:n] now contains sum of hits for each seq, but we want just the increment
        // lcounts[0:n] now contains count of 1.0 hits for each seq, but we want just the increment
        int start=0;
        for (int specc=firstSpectraInChunk; specc<firstSpectraNotInChunk; specc++) {
            int end = sequenceCountAfterEachScorePreload[specc] - devOffset;
            thrust::adjacent_difference(
                miUsed.begin()+start,
                miUsed.begin()+end,
                dScoresDev.begin()+start+devOffset);
            thrust::adjacent_difference(
                lcounts.begin()+start,
			    lcounts.begin()+end,
			    lCountsDev.begin()+start+devOffset);
            start = end;
        }
        devOffset += start; // in case we're looping
#ifdef ___DEBUG
STDVECT(int,scatteredCopiesLocal,lCountsDev);
printf("lcounts ");for (int nnn=0;nnn<scatteredCopiesLocal.size();nnn++)printf("%d,%d ",nnn,scatteredCopiesLocal[nnn]);printf("\n");fflush(stdout);
STDVECT(float,dscoreLocal,dScoresDev);
printf("dscores ");for (int nnn=0;nnn<dscoreLocal.size();nnn++)printf("%d,%f ",nnn,dscoreLocal[nnn]);printf("\n\n");fflush(stdout);
int blah=0;blah++;
#endif
    }
    // copy back to host memory
    mscore_kgpu_thrust_device_to_host_copy_int(lCountsDev,lCountsResult);

    // copy back to host memory
    mscore_kgpu_thrust_device_to_host_copy_float(dScoresDev,dScoresResult);

    cudaStreamSynchronize(MY_CUDA_STREAM); // wait for memcpy to complete

    return;
}

