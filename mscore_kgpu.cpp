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

mscore_kgpu::mscore_kgpu(void) : mscore_k(), m_miUsed(NULL)
{
}

mscore_kgpu::~mscore_kgpu(void)
{
    clear();
    mscore_kgpu_thrust_fvec_kill(m_miUsed);
}

bool mscore_kgpu::clear()
{
    mscore_k::clear();
    for (int n=m_vmiTypeGPU.size();n--;) {
        m_vmiTypeGPU[n].kill();
    }
    m_vmiTypeGPU.resize(0);
    return true;
}

/*
 * called before scoring inside the score() function to allow any
 * necessary resetting of member variables.
 */
void mscore_kgpu::prescore(const size_t _i)
{
    mscore::prescore(_i);

    // Initialize of clear the used intensity look-up.
    if (m_miUsed == NULL)
        m_miUsed = mscore_kgpu_thrust_fvec_alloc(m_maxEnd);
    else
        mscore_kgpu_thrust_fvec_clear(m_miUsed);
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

    vmiTypeGPU vTypeGPU;
    vTypeGPU.init(_s.m_vMI.size());
    if (_s.m_vMI.size() != 0)
    {
        // CUDA works better on arrays than structs, lay these out
        thrust::host_vector<float> hfI; // intensities
        thrust::host_vector<float> hfM; // m/z's
        hfI.resize( _s.m_vMI.size());
        hfM.resize( _s.m_vMI.size());
        for (int n=_s.m_vMI.size();n--;) {
            hfI[n] = _s.m_vMI[n].m_fI;
            hfM[n] = _s.m_vMI[n].m_fM;
        }
        // Screen peeks on upper end.
        int iWindowCount = 10;
        int endMassMax = (int)(((_s.m_dMH + (_s.m_fZ - 1) * m_seqUtil.m_dProton) / _s.m_fZ) * 2.0 + 0.5) + iWindowCount;

        // now pass off to CUDA implementation
        mscore_kgpu_thrust_score(hfI,hfM,m_dIsotopeCorrection,iWindowCount,endMassMax,m_maxEnd,vTypeGPU); // results come back in vTypeGPU
    }
    m_vmiTypeGPU.push_back(vTypeGPU);
    return true;

}


/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses).
 */
double mscore_kgpu::dot(unsigned long *_v)
{
    return mscore_kgpu_thrust_dot(*_v,m_vmiTypeGPU[m_lId],m_plSeq,*m_miUsed);
}

