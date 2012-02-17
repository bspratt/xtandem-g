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

#include "mscore_kgpu.h"
#include "cuda.h"
#include <iostream>
#include <sstream>

// using CUDA Thrust template lib
#include "mscore_kgpu_thrust.h"
#include <thrust/experimental/cuda/pinned_allocator.h> 

// Factory instance, registers itself with the mscoremanager.
static mscorefactory_kgpu factory;
static bool first_access = true;



// format a byte count in a human readable manner
std::string prettymem(size_t bytes) {
  std::ostringstream result;
  result << bytes << " bytes";
  if (bytes >= (1024*1024)) {
	result << " (" << bytes/(1024*1024) <<" MB)";
  }
  else if (bytes >= 1024 ) {
	result << " (" <<bytes/1024 <<" KB)";
  }
  return result.str();
}

// logging stuff, stolen from FastPass project
enum LogType {logall, debug, info, warn, error, lognone};
enum EventType {noevent, init, initok, hostalloc, hostrelease, hostreleaseok, hostallocok, devicealloc, devicerelease, devicereleaseok, deviceallocok, parse, parseok, copye, copyok, kernel, kernelok, file, fileok};

// ug, better design is wanted here
static LogType logLevel=debug;

#define logevent(LOGTYPE, EVENT, LOGOUTPUT)	\
  if (LOGTYPE >= logLevel) { \
  std::cerr << #LOGTYPE << "(" << #EVENT << "): " << LOGOUTPUT ; \
  }

// end logging stuff


mscorefactory_kgpu::mscorefactory_kgpu()
{
    m_initialFreeMemory = 0;
    mscoremanager::register_factory("kgpu-score", this);
}

mplugin* mscorefactory_kgpu::create_plugin()
{
    if (first_access) { // first call, from master thread
        first_access = false; // just once, thanks
        // initialize GPU
        int deviceCount = 0;
        int bestDeviceIndex = 0;
        double highestDeviceScore = 0.0;  // device score = clockRate * multiProcessorCount
        cudaDeviceProp deviceProp;

        cuInit(0);

        // loop through all devices and choose which one to use based on
        // highest device score

        cudaError_t err = cudaGetDeviceCount(&deviceCount);
        if (err != cudaSuccess) {
             logevent(error, init, cudaGetErrorString(err) << endl);
        }
        logevent(info, init, "GPU devices found: " << deviceCount << endl);

        for (int curDeviceNumber=0; curDeviceNumber < deviceCount; curDeviceNumber++) {
            logevent(info, init, "testing device " << curDeviceNumber << endl);
            memset(&deviceProp, 0, sizeof(deviceProp));
            if (cudaSuccess == cudaGetDeviceProperties(&deviceProp, curDeviceNumber)) {
                double curDeviceScore = deviceProp.multiProcessorCount * deviceProp.clockRate;
                if (curDeviceScore > highestDeviceScore) {
                    highestDeviceScore = curDeviceScore;
                    bestDeviceIndex = curDeviceNumber;
                }
            }
            else {
                logevent(error, init, "unable to access GPU properites, using non-GPU k-score" << endl);
                return new mscore_k();
            }
        }

        memset(&deviceProp, 0, sizeof(deviceProp));
        if (cudaSuccess == cudaGetDeviceProperties(&deviceProp, bestDeviceIndex)) {
            logevent(info, init, "Using GPU device " << bestDeviceIndex << ": \"" << deviceProp.name << "\", "<< prettymem(deviceProp.totalGlobalMem) << ", " << deviceProp.multiProcessorCount << " multiprocessors, " << deviceProp.clockRate/1000 << " MHz" << endl);
            size_t free_mem, total_mem;
            cudaMemGetInfo(&free_mem, &total_mem);
            m_initialFreeMemory = free_mem;
            logevent(info, init, "device initial free memory: " << prettymem(m_initialFreeMemory) << endl);
        }
        else {
            logevent(error, init, "unable to access device properites, using non-GPU k-score" << endl);
            return new mscore_k();
        }      

        logevent(info, init, "resetting GPU" << endl);
        if (cudaSuccess != cudaDeviceReset() ) {
            logevent(error, init, "error resetting GPU, using non-GPU k-score" << endl);
            return new mscore_k();
        }
        logevent(info, initok, "done resetting GPU" << endl);

        if (cudaSuccess != cudaSetDevice(bestDeviceIndex)) {
            logevent(error, init, "Error - cannot set GPU, using non-GPU k-score"  << endl);
            return new mscore_k();
        }
        else {
            logevent(info, initok, "GPU " << bestDeviceIndex << " initialized." << endl);
        }
        size_t free, total;
        free = total = 0;
        cudaMemGetInfo(&free, &total);
        m_initialFreeMemory = total;
        logevent(info, init, "Device memory: " << prettymem(free) << " free / " << prettymem(m_initialFreeMemory) << " total" << endl);

    }
    if (!m_initialFreeMemory) {
        logevent(error, init, "Error - insufficient GPU memory, using non-GPU k-score" << endl);
        return new mscore_k();
    }
    mscore_kgpu_thrust_init();
    return new mscore_kgpu();
}

mscore_kgpu::mscore_kgpu(void) : mscore_k(), m_preloading(false), m_cached_sequences_i(NULL)
{
    clear();
}

mscore_kgpu::~mscore_kgpu(void)
{
    clear();
  
}

bool mscore_kgpu::clear()
{
    bool ret = mscore_k::clear();
    for (int n=m_vSpectraGPU.size();n--;) {
        m_vSpectraGPU[n].kill();
    }
    m_vSpectraGPU.resize(0);
    clear_sequence_cache();
    return ret;
}

void mscore_kgpu::clear_sequence_cache() {
    m_cached_sequences_l.resize(0);
    m_cached_sequences_index.resize(0);
    m_cached_sequence_lengths.resize(0);
    m_previous_lId = -1;
}

bool mscore_kgpu::load_next(void) {  
    return mscore_k::load_next();  // TODO - lookahead so we can gang up dot() work

#ifdef TODO
     needs to preprocess this logic in the calling mprocess object and play back the m_pScore behavior:

    			    while(m_pScore->load_next())	{
						m_tPeptideCount += m_pScore->m_State.m_lEqualsS;
						m_bPermute = false;
						m_bPermuteHigh = false;
						create_score(_s,lStart,lEnd,lMissedCleaves - 1,true);
						if(m_bPermute && m_bCrcCheck || m_bPermuteHigh)	{
							m_pScore->reset_permute();
							while(m_pScore->permute())	{
								create_score(_s,lStart,lEnd,lMissedCleaves - 1,false);
							}
						}
					}

which means reproducing the minumum functionality of mprocess::create_score and playing back the m_pScore behavior

bool mprocess::create_score(const msequence &_s,const size_t _v,const size_t _w,const long _m,bool _p)
{
	long lIonCount = 0;
	float fScore = -1.0;
	float fHyper = -1.0;
	size_t a = 0;
	size_t b = 0;
	size_t c = 0;
	bool bOk = false;
	bool bDom = false;
	long lCount = 0;
	bool bIonCheck = false;
	bool bMassCheck = false;
/*
 * score each mspectrum identified as a candidate in the m_pScore->m_State object
 */
	while(lCount < m_pScore->m_State.m_lEqualsS)	{
		a = m_pScore->m_State.m_plEqualsS[lCount];
		lCount++;
		lIonCount = 0;
/*
* this check is needed to keep tandem consistent whether running on
* single-threaded, multi-threaded or on a cluster.  otherwise, when
* there are multiple spectra matching a sequence (which is more common
* the fewer mprocess objects there are) one can cause others to score
* more permutation sequences.
*/
		if (!_p && m_vSpectra[a].m_hHyper.m_ulCount >= 400)
			continue;
		fScore = 1.0;
		fHyper = 1.0;
		m_pScore->m_lMaxCharge = (long)(m_vSpectra[a].m_fZ+0.1);
/*
 * in versions prior to 2004.03.01, spectra with m_bActive == false were
 * rejected at this point, to save time & because of a problem with
 * multiple recording of the same sequence. starting with 2004.03.01,
 * the later problem has been corrected, and because of point mutation
 * analysis, it has become important to reexamine all sequences.
 */
		fScore = m_pScore->score(a);
		fHyper = m_pScore->m_fHyper;
    }
#endif
}


/*
 * called before scoring inside the score() function to allow any
 * necessary resetting of member variables.
 */
void mscore_kgpu::prescore(const size_t _i)
{
    mscore_k::prescore(_i);
    if (m_preloading) {  // avoid infinite recursion
        return;
    }

    // check to see if we can just use the same sequence modifications again
   m_sequenceIsSameAsLastTime = false;
    if ((m_previous_lId>=0) && (this->m_vSpec[m_previous_lId].m_fZ == m_vSpec[m_lId].m_fZ) &&
        (m_previousSequence==this->m_pSeq)) {
        int n;
        long maxSeqMZ = (long)m_vSpectraGPU[m_lId].iM_max; // ignore any sequence members greater that largest spectra mz
        for (n=m_cached_sequences_index.size();n-->1;) {
            if (m_cached_sequences_l[m_cached_sequences_index[n]-1] > maxSeqMZ) {
                break;
            }
        }
        if (!n) {
            m_sequenceIsSameAsLastTime = true;
        }
    }
    if (m_sequenceIsSameAsLastTime) {
        // we can just play back the same sequence modifications
        m_current_playback_sequence = 0;
    } else {
        // preprocess data for dot() - runs
        // mscore::score, without the calls to dot()
        m_previousSequence = m_pSeq;
        n_cached_sequences = 0;
        clear_sequence_cache();
        m_vSpectra.clear();
        m_preloading = true;
        score(_i); // runs through gathering sequences without performing dot()
        m_preloading = false;
        m_current_playback_sequence = 0;
        //  convert to int and push to device memory
        pinned_host_vector_int_t cached_sequences_h(m_cached_sequences_l.size());
        pinned_host_vector_int_t::iterator cs = cached_sequences_h.begin();
        for (std::vector<long>::iterator it = m_cached_sequences_l.begin(); it != m_cached_sequences_l.end();) {
            *cs++ = *it++;
        }
        m_cached_sequences_i = mscore_kgpu_thrust_device_copy(cached_sequences_h,m_cached_sequences_i);
    }
    m_previous_lId = m_lId; // next time we come through here, see if we can reuse all this
}

bool mscore_kgpu::load_seq(const unsigned long _t,const long _c) {
    if (m_preloading) {
        bool ret = mscore_k::load_seq(_t,_c); // get it
        cache_sequence();
        return ret;
    } else {
        return playback_sequence();
    }
}

void mscore_kgpu::cache_sequence() {
    m_cached_sequence_lengths.push_back(m_lSeqLength); // preserve state - not easily derivable

    // jam sequence info into contiguous memory for cheap transfer to device
    int seqlen;
    // get length of new sequence, ignore any members greater than max spectrum mz
    for (seqlen=0 ; m_plSeq[seqlen] && (m_plSeq[seqlen] <= (unsigned long)m_vSpectraGPU[m_lId].iM_max); seqlen++); 
    int prev_end=m_cached_sequences_l.size(); // note end of list
    m_cached_sequences_index.push_back(prev_end); // note where to find new sequence
    m_cached_sequences_l.resize(m_cached_sequences_l.size()+seqlen);
    memmove(&m_cached_sequences_l[prev_end],m_plSeq,seqlen*sizeof(long));
    n_cached_sequences++;
}

bool mscore_kgpu::playback_sequence() {
    m_lSeqLength = m_cached_sequence_lengths[m_current_playback_sequence];
    m_current_playback_sequence_begin = m_cached_sequences_index[m_current_playback_sequence++];
    m_current_playback_sequence_end = (m_current_playback_sequence==m_cached_sequences_index.size())?
        m_cached_sequences_l.size():
        m_cached_sequences_index[m_current_playback_sequence];
   int len = (m_current_playback_sequence_end-m_current_playback_sequence_begin);
    memmove(m_plSeq,&m_cached_sequences_l[m_current_playback_sequence_begin], len*sizeof(long));
    m_plSeq[len] = 0;
    return true;
}

bool mscore_kgpu::add_A(const unsigned long _t,const long _c) {
    if (m_preloading) {
        return mscore_k::add_A(_t,_c); // get it
    }
    return true;
}
bool mscore_kgpu::add_B(const unsigned long _t,const long _c) {
    if (m_preloading) {
        return mscore_k::add_B(_t,_c); // get it
    }
    return true;
}
bool mscore_kgpu::add_C(const unsigned long _t,const long _c) {
    if (m_preloading) {
        return mscore_k::add_C(_t,_c); // get it
    }
    return true;
}
bool mscore_kgpu::add_Y(const unsigned long _t,const long _c) {
    if (m_preloading) {
        bool ret = mscore_k::add_Y(_t,_c); // get it
        if (!_t) { // direct call   
            cache_sequence();
        }
        return ret;
    } else {
         if (!_t) { // direct call
            playback_sequence();
         }
    }
    return true;
}
bool mscore_kgpu::add_X(const unsigned long _t,const long _c) {
    if (m_preloading) {
        return mscore_k::add_X(_t,_c); // get it
    }
    return true;
}
bool mscore_kgpu::add_Z(const unsigned long _t,const long _c) {
    if (m_preloading) {
        return mscore_k::add_Z(_t,_c); // get it
    }
    return true;
}


/*
 * add_mi does the work necessary to set up an mspectrum object for modeling. 
 *   - an entry in the m_State object is made for the parent ion M+H
 * once an mspectrum has been added, the original mspectrum is no longer
 * needed for modeling, as all of the work associated with a spectrum
 * is only done once, prior to modeling sequences.
 */
bool mscore_kgpu::add_mi(mspectrum &_s)
{
    if (!mscore::add_mi(_s))
        return false;

    if (&_s == m_vSpectra.back()) { // last in list, transfer all to device memory in one go
cudatimer mi_time;CUDA_TIMER_START(mi_time);
        size_t a, n_pairs=0;
        for (a=0;a < m_vSpectra.size();a++)	{
            n_pairs += m_vSpectra[a]->m_vMI.size();
        }
        pinned_host_vector_float_t fM(n_pairs);
        pinned_host_vector_float_t fI(n_pairs);
        n_pairs = 0;
        for (a=0;a < m_vSpectra.size();a++)	{
            for (size_t n=0;n<m_vSpectra[a]->m_vMI.size();n++) {
                fI[n_pairs]   = m_vSpectra[a]->m_vMI[n].m_fI;
                fM[n_pairs++] = m_vSpectra[a]->m_vMI[n].m_fM;
            }
        }
        // now copy to memory
        thrust::device_vector<float> dfM(fM);
        thrust::device_vector<float> dfI(fI);

        // and actually do the add_mi logic
        n_pairs = 0;
        for (a=0;a < m_vSpectra.size();a++)	{
            mspectrum &_s = *m_vSpectra[a];
            vmiTypeGPU vTypeGPU;
            vTypeGPU.init(_s.m_vMI.size());
            if (_s.m_vMI.size() != 0)
            {
               // Screen peeks on upper end.
                int iWindowCount = 10;
                int endMassMax = (int)(((_s.m_dMH + (_s.m_fZ - 1) * m_seqUtil.m_dProton) / _s.m_fZ) * 2.0 + 0.5) + iWindowCount;

                // now pass off to CUDA implementation
                mscore_kgpu_thrust_score(dfM.begin()+n_pairs,dfI.begin()+n_pairs,_s.m_vMI.size(),m_dIsotopeCorrection,iWindowCount,endMassMax,m_maxEnd,vTypeGPU); // results come back in vTypeGPU
            }
            m_vSpectraGPU.push_back(vTypeGPU);
            n_pairs += _s.m_vMI.size();
        }
CUDA_TIMER_STOP(mi_time)
    }
    return true;

}

bool mscore_kgpu::add_details(mspectrum &_s)
{
    bool ret = mscore::add_details(_s); // do the base work
    m_vSpectra.push_back(&_s); // note who we'll be working on
    return ret;
}

/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses).
 */
int kgfoo = 0;
double mscore_kgpu::dot(unsigned long *_v)
{
    if (m_preloading) {
        return 0.0;  // not yet, still loading spectra
    } else {
        if (!m_current_playback_sequence_begin) {
            if (!m_sequenceIsSameAsLastTime) {
                m_dScores.resize(n_cached_sequences);
                m_lCounts.resize(n_cached_sequences);
                m_cached_sequences_index.push_back(m_cached_sequences_l.size()); // note end
            }
            mscore_kgpu_thrust_dot(m_lCounts,m_dScores,m_vSpectraGPU[m_lId],
                m_cached_sequences_i->begin(),m_cached_sequences_index,
                m_sequenceIsSameAsLastTime);
        }


        *_v = (unsigned long)m_lCounts[m_current_playback_sequence-1];
if (!(kgfoo++%1000)) std::cout<<kgfoo<<"\n";
//printf("%d %d %f\n",kgfoo,*_v,m_dScores[m_current_playback_sequence-1]);
        return m_dScores[m_current_playback_sequence-1];
    }
}

