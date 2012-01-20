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

class mscore_kgpu : public mscore_k
{
#ifdef HAVE_MULTINODE_TANDEM // support for Hadoop and/or MPI?
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<mscore>(*this);
      ar & m_maxEnd;
      ar & m_miUsed;
      ar & m_vmiType;
      ar & m_dIsotopeCorrection;
    }
#endif // HAVE_MULTINODE_TANDEM
protected:
    friend class mscorefactory_kgpu;

    mscore_kgpu(void);    // Should only be created through mscorefactory_tandem

public:
    virtual ~mscore_kgpu(void);

public:
    virtual bool load_param(XmlParameter &_x); // allows score object to issue warnings,
                                               // or set variables based on xml
    virtual bool precondition(mspectrum &_s); // called before spectrum conditioning
    virtual void prescore(const size_t _i); // called before scoring
    
    virtual bool add_mi(mspectrum &_s);

    virtual double sfactor(); // factor applied to final convolution score
    virtual unsigned long mconvert(double _m, const long _c); // convert mass to integer ion m/z for mi vector
    virtual void report_score(char* _buff, float _h); // format hyper score for output

    virtual bool clear();

protected:
    virtual double dot(unsigned long *_v); // this is where the real scoring happens

protected:
    unsigned long imass(double _m)
    {
        return (unsigned long)((_m/m_dIsotopeCorrection) + 0.5);
    }

protected:
    int m_maxEnd;
    miLookup m_miUsed;
    vectorvmiType m_vmiType;
		vector<int> m_vmiUsed;

    double m_dIsotopeCorrection;
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
