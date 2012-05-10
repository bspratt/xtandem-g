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

#ifndef MSCORE_KGPU_THRUST_H
#define MSCORE_KGPU_THRUST_H


// using CUDA Thrust template lib
#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/experimental/cuda/pinned_allocator.h> 

#include "mspectrum.h"

typedef thrust::host_vector<
    float,
    thrust::experimental::cuda::pinned_allocator<float> 
  > pinned_host_vector_float_t;

typedef thrust::host_vector<
    int,
    thrust::experimental::cuda::pinned_allocator<int> 
  > pinned_host_vector_int_t;


class cudatimer {
 public: 
  cudaEvent_t start;
  cudaEvent_t stop;
  float elapsed;
  float count;
  float period;
 cudatimer() : elapsed(0),count(0),period(1000) { }
};

#if 1 || defined(PROFILE)
#define CUDA_TIMER_DECL(tvar)
#define CUDA_TIMER_START(tvar)
#define CUDA_TIMER_STOP(tvar) 
#define myassert(expr) 
#else
#define CUDA_TIMER_DECL(tvar) cudatimer tvar;
#define CUDA_TIMER_START(tvar)      \
  cudaEventCreate(&tvar.start);	    \
  cudaEventCreate(&tvar.stop);      \
  cudaEventRecord( tvar.start, 0 );

#define CUDA_TIMER_STOP(tvar) \
  cudaEventRecord( tvar.stop, 0 ); \
  cudaEventSynchronize( tvar.stop ); \
  {float t;cudaEventElapsedTime( &t, tvar.start, tvar.stop ); tvar.elapsed+=t;} \
  if (!(((int)++tvar.count)%((int)tvar.period))) \
    std::cout << "t("<<#tvar<<")="<<tvar.elapsed/(tvar.count)<<"ms avg\n"; \
  cudaEventDestroy( tvar.start ); cudaEventDestroy( tvar.stop ); \

#define myassert(expr) \
  if (!(expr))								\
    {std::cerr<<"assert failed at "<< __FILE__<<":"<<__LINE__<<", " \
	      << #expr <<"\n";exit(1);}
#endif


// helpful macro for copying to c++ space for debugger viewing
#ifdef _DEBUG
#define STDVECT(T,local,mem)  \
  std::vector<T> local; \
  for (size_t nn=0;nn<mem.size();)local.push_back(mem[nn++]);
#else
#define STDVECT(T,local,mem) 
#endif

void mscore_kgpu_thrust_init(); // call once at start

// helpers for hiding device memory implementation from normal C++
thrust::device_vector<float> *mscore_kgpu_thrust_fvec_alloc(int size);
void mscore_kgpu_thrust_fvec_clear(thrust::device_vector<float> *vec);
void mscore_kgpu_thrust_fvec_kill(thrust::device_vector<float> *vec);
void mscore_kgpu_thrust_ivec_kill(thrust::device_vector<int> *vec);
thrust::device_vector<int> *mscore_kgpu_thrust_host_to_device_copy(
			      const pinned_host_vector_int_t &host_src,
			      thrust::device_vector<int> *device_dest);
void mscore_kgpu_thrust_host_to_device_copy_float(
  const pinned_host_vector_float_t &host_src,
  thrust::device_vector<float> &device_dest);
void mscore_kgpu_thrust_device_to_host_copy_float(
  const thrust::device_vector<float> &device_src,
  pinned_host_vector_float_t &host_dest);
void mscore_kgpu_thrust_device_to_host_copy_int(
  const thrust::device_vector<int> &device_src,
  pinned_host_vector_int_t &hhost_dest);


struct vmiTypeGPU {
  // purposefully crude, don't want lots of automatic create destroy
  // copy in std::vector
public:
    thrust::device_vector<int> *iM; // m/z's converted to integer indexes
    int iM_max; // iM->back() in c++ space for cheap access
    thrust::device_vector<float> *fI; // intensities: for ion=iM[n],
				      // intensty=fI[n]
    void init(int s);
    void kill();
};

void mscore_kgpu_thrust_score(
    // m/z-intensity pairs (input)
    const thrust::device_vector<float>::iterator fM,
    const thrust::device_vector<float>::iterator fI,    
    size_t length, // iterator end distance
    double dIsotopeCorrection, // for m/z binning (input)
    int iWindowCount, // (input)
    int endMassMax, // max acceptable binned m/z (input)
    int &m_maxEnd,  // max encountered binned m/z (output)
    vmiTypeGPU &binned);   // binned intensities and m/z's (output)


/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses).
 */
// perform dot on current spectrum and all sequence variations at one go
void mscore_kgpu_thrust_dot(pinned_host_vector_int_t &lCountsResult,
    pinned_host_vector_float_t &dScoresResult,
    const std::vector<const vmiTypeGPU *> &spectra,
    const std::vector<int> &sequenceCountAfterEachScorePreload,
    const thrust::device_vector<int> &cached_sequences, 
    const std::vector<int> &sequence_index);

#endif // MSCORE_KGPU_THRUST_H

