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
    return new mscore_kgpu();
}

mscore_kgpu::mscore_kgpu(void) : mscore_k(), m_miUsed(NULL), m_preloading(false), m_cached_sequences_i(NULL)
{
}

mscore_kgpu::~mscore_kgpu(void)
{
    clear();
    mscore_kgpu_thrust_fvec_kill(m_miUsed);
    
}

bool mscore_kgpu::clear()
{
    bool ret = mscore_k::clear();
    for (int n=m_vmiTypeGPU.size();n--;) {
        m_vmiTypeGPU[n].kill();
    }
    m_vmiTypeGPU.resize(0);
    clear_sequence_cache();
    return ret;
}

void mscore_kgpu::clear_sequence_cache() {
    mscore_kgpu_thrust_ivec_kill(m_cached_sequences_i);
    m_cached_sequences_i = NULL;
    m_cached_sequences_l.resize(0);
    m_cached_sequences_f.resize(0);
    m_cached_sequences_index.resize(0);
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

    // Initialize of clear the used intensity look-up.
    if (m_miUsed == NULL)
        m_miUsed = mscore_kgpu_thrust_fvec_alloc(m_maxEnd);
    else
        mscore_kgpu_thrust_fvec_clear(m_miUsed);
    clear_sequence_cache();

    // preprocess data for dot() - runs
    // mscore::score, without the calls to dot()
    m_vSpectra.clear();
    m_preloading = true;
    score(_i); // runs through gathering sequences without performing dot()
    m_preloading = false;
    m_current_playback_sequence = 0;
    //  convert to in and push to device memory
    pinned_host_vector_int_t cached_sequences_h(m_cached_sequences_l.size());
    pinned_host_vector_int_t::iterator cs = cached_sequences_h.begin();
    for (std::vector<long>::iterator it = m_cached_sequences_l.begin(); it != m_cached_sequences_l.end();) {
        *cs++ = *it++;
    }
    m_cached_sequences_i = mscore_kgpu_thrust_device_copy(cached_sequences_h);
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
    // jam sequence info into contiguous memory for cheap transfer to device
    int seqlen;
    for (seqlen=0 ; m_plSeq[seqlen]; seqlen++); // get length of new sequence
    int prev_end=m_cached_sequences_l.size(); // note end of list
    m_cached_sequences_index.push_back(prev_end); // note where to find new sequence
    m_cached_sequences_f.resize(m_cached_sequences_f.size()+seqlen+1);
    memmove(&m_cached_sequences_f[prev_end],this->m_pfSeq,(seqlen+1)*sizeof(float));
    m_cached_sequences_l.resize(m_cached_sequences_l.size()+seqlen+1);
    memmove(&m_cached_sequences_l[prev_end],this->m_plSeq,(seqlen+1)*sizeof(long));
}

bool mscore_kgpu::playback_sequence() {
    m_current_playback_sequence_begin = m_cached_sequences_index[m_current_playback_sequence++];
    m_current_playback_sequence_end = (m_current_playback_sequence==m_cached_sequences_index.size())?
        m_cached_sequences_l.size():
        m_cached_sequences_index[m_current_playback_sequence];
    memmove(m_pfSeq,&m_cached_sequences_f[m_current_playback_sequence_begin],(m_current_playback_sequence_end-m_current_playback_sequence_begin)*sizeof(float));
    memmove(m_plSeq,&m_cached_sequences_l[m_current_playback_sequence_begin],(m_current_playback_sequence_end-m_current_playback_sequence_begin)*sizeof(long));
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
                mscore_kgpu_thrust_score(dfM,dfI,n_pairs,n_pairs+_s.m_vMI.size(),m_dIsotopeCorrection,iWindowCount,endMassMax,m_maxEnd,vTypeGPU); // results come back in vTypeGPU
            }
            m_vmiTypeGPU.push_back(vTypeGPU);
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
double mscore_kgpu::dot(unsigned long *_v)
{
    if (m_preloading) {
        return 1.0;  // not yet, still loading spectra
    } else {
        return mscore_kgpu_thrust_dot(*_v,m_vmiTypeGPU[m_lId],m_cached_sequences_i, m_current_playback_sequence_begin, m_current_playback_sequence_end,*m_miUsed);
    }
}

