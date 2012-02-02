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
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

// helpers for hiding device memory implementation from normal C++
thrust::device_vector<float> *mscore_kgpu_thrust_fvec_alloc(int size);
void mscore_kgpu_thrust_fvec_clear(thrust::device_vector<float> *vec);
void mscore_kgpu_thrust_fvec_kill(thrust::device_vector<float> *vec);

struct vmiTypeGPU {
    // purposefully crude, don't want lots of automatic create destroy copy in std::vector
public:
    thrust::device_vector<float> *fI; // m/z's converted to integer indexes
    thrust::device_vector<int> *iM; // m/z's converted to integer indexes
    void init(int s);
    void kill();
};

void mscore_kgpu_thrust_score(thrust::host_vector<float> const &host_fI, // intensities (input)
    thrust::host_vector<float> const &host_fM, // m/z's (input)
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
float mscore_kgpu_thrust_dot(unsigned long &lCount,const vmiTypeGPU &_vmiTypeGPU,
    unsigned long const *_plSeq, thrust::device_vector<float> &_miUsed);

#endif // MSCORE_KGPU_THRUST_H

