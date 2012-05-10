// -*- mode: c++ -*-

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
  void serialize(Archive &ar, const unsigned int version) {
    ar & boost::serialization::base_object<mscore_k>(*this);
  }
#endif // HAVE_MULTINODE_TANDEM
protected:
    friend class mscorefactory_kgpu;

    mscore_kgpu(void);    // Should only be created through mscorefactory_tandem

    pinned_host_vector_float_t m_dScores;
    pinned_host_vector_int_t m_lCounts;
    std::vector<vmiTypeGPU> m_vSpectraGPU; // processed intensity-m/z pairs
    std::vector<mspectrum *> m_vSpectraToScore; // spectra we'll be scoring

public:
    virtual ~mscore_kgpu(void);

public: // inherit anything we don't cuda-ize
    virtual bool load_next(const mprocess *_parentProcess); 
    virtual void prescore(const size_t _i);
    virtual double dot(unsigned long *_v); // this is where the real
					   // scoring happens
    virtual bool add_mi(mspectrum &_s);
    virtual bool add_details(mspectrum &_s);
    virtual bool clear();

    // stuff for batching up sequence transfers to device memory
    bool m_preloading;
    virtual bool load_seq(const unsigned long _t,const long _c);
protected:
    virtual bool add_A(const unsigned long _t,const long _c);
    virtual bool add_B(const unsigned long _t,const long _c);
    virtual bool add_C(const unsigned long _t,const long _c);
    virtual bool add_Y(const unsigned long _t,const long _c);
    virtual bool add_X(const unsigned long _t,const long _c);
    virtual bool add_Z(const unsigned long _t,const long _c);

    // for lookahead in mscore::load_next
    bool m_FakeProcess_bPermute;
    bool m_FakeProcess_bPermuteHigh;
    bool fakeProcess_create_score(const mprocess *_parentProcess,bool _p);
    inline void assert_consistency(); // no effect unless _DEBUG is defined
    void clear_sequence_cache();
    void cache_sequence();
    bool playback_sequence();
    std::vector<int> m_sequenceCountAfterEachScorePreload;
    std::vector<const vmiTypeGPU *> m_currentSpectra;
    int m_current_playback_sequence;
    int m_current_playback_sequence_begin;
    int m_current_playback_sequence_end;
    pinned_host_vector_int_t m_cached_sequences_l;
    thrust::device_vector<int> *m_cached_sequences_i;

    class mscore_internals_cacheinfo {
    public:
        mscore_internals_cacheinfo(){};
        mscore_internals_cacheinfo(int sequence_length, double seqMH, int lId,  
            float fMinMass, float fMaxMass,
            unsigned long lMaxPeaks, long lMaxCharge,
            const char *seq_p, float *fseq_p,unsigned long *lseq_p, 
            const unsigned long *plCount,
        	const float *pfScore, // convolute score information,
				      // indexed using the
				      // mscore_type_a enum
            unsigned long lType);
        int m_lSeqLength;
        double m_seqMH; 
        int m_lId;
        float m_fMinMass;
        float m_fMaxMass;
        long m_lMaxCharge; // current parent ion charge
        unsigned long m_lMaxPeaks; // if > 0, the m_lMaxPeaks most intense peaks will be used
        std::string m_sequence;
        std::vector<float> m_fseq;
        std::vector<long> m_lseq;
	unsigned long m_plCount[16]; // ion count information, indexed
				     // using the mscore_type_a enum
	float m_pfScore[16]; // convolute score information, indexed
			     // using the mscore_type_a enum
        unsigned long m_lType; // current ion type - value from mscore_type
    };
    std::vector<mscore_internals_cacheinfo *> m_mscore_internals_cache;
    void save_mscore_internals() {
        m_mscore_internals_cache.
	  push_back(new mscore_internals_cacheinfo(
		      m_lSeqLength, m_dSeqMH,m_lId, m_fMinMass, m_fMaxMass, 
		      m_lMaxPeaks, m_lMaxCharge, m_pSeq, m_pfSeq, m_plSeq, 
		      m_plCount, m_pfScore, m_lType));
    }
    void restore_mscore_internals(int n) {
        const mscore_internals_cacheinfo &c=*m_mscore_internals_cache[n];
        m_lSeqLength = c.m_lSeqLength;
        m_lId = c.m_lId;
	    m_fMinMass = c.m_fMinMass;
	    m_fMaxMass = c.m_fMaxMass;
        m_lMaxPeaks = c.m_lMaxPeaks; 
        m_lMaxCharge = c.m_lMaxCharge;
        m_lType = c.m_lType;
        if (m_pSeq) {
            strcpy(m_pSeq,c.m_sequence.c_str());
        }
        if (m_pfSeq) {
            memmove(m_pfSeq,&c.m_fseq[0],sizeof(float)*c.m_fseq.size());
        }
        if (m_plSeq) {
            memmove(m_plSeq,&c.m_lseq[0],sizeof(long)*c.m_lseq.size());
        }
        m_dSeqMH = c.m_seqMH;
        memmove(m_plCount,c.m_plCount,16*sizeof(long));
        memmove(m_pfScore,c.m_pfScore,16*sizeof(float));
    }
    std::vector<int> m_cached_sequences_index;
    std::vector<mscorestate> m_cached_mscorestates;
    int m_currentStateIndex;
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
