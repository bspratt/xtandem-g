/*
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
/*
  modified December 2010 by Insilicos LLC to support object serialization for 
  Hadoop and MPI use
*/
#ifndef MSCORE_KGPU_H
#define MSCORE_KGPU_H

#include <cuda_runtime.h> // must include before stdafx.h to avoid redefines
#include "stdafx.h"
#include "msequence.h"
#include "mspectrum.h"
#include "msequtilities.h"
#include "xmlparameter.h"
#include "mscore_k.h" // we want to be like k-score, but faster
#include "mscore_kgpu_thrust.h" // stuff that gets compiled by NVCC

class mscore_kgpu : public mscore_k
{
#ifdef HAVE_MULTINODE_TANDEM // support for Hadoop and/or MPI?
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<mscore_k>(*this);
    }
#endif // HAVE_MULTINODE_TANDEM
protected:
    friend class mscorefactory_kgpu;

    mscore_kgpu(void);    // Should only be created through mscorefactory_tandem

    thrust::device_vector<float> *m_miUsed; // intensities of used masses
    std::vector<vmiTypeGPU> m_vmiTypeGPU; // processed intensity-m/z pairs

public:
    virtual ~mscore_kgpu(void);

public: // inherit anything we don't cuda-ize
    virtual void prescore(const size_t _i);
    virtual double dot(unsigned long *_v); // this is where the real scoring happens
    virtual bool add_mi(mspectrum &_s);
    virtual bool clear();
};

/*
 * mscorefactory_kgpu implements a factory for creating mscore_k instances.
 */
class mscorefactory_kgpu : public mpluginfactory
{
public:
    mscorefactory_kgpu();
    // GPU stats and info
    size_t m_initialFreeMemory;

    virtual mplugin* create_plugin();
};

#endif // MSCORE_KGPU_H
