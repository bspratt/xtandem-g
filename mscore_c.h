/* -*-c++-*-
 *
 * Copyright (c) 2008-2010 Fred Hutchinson Cancer Research Center
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
#ifndef MSCORE_C_H
#define MSCORE_C_H

#include <cassert>
#include "mscore.h"

typedef vector<mi>::const_iterator vmiIter;

class mscore_c : public mscore
{
#ifdef HAVE_MULTINODE_TANDEM // support for Hadoop and/or MPI?
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
	  ar & boost::serialization::base_object<mscore>(*this);
	  ar & m_miUsed;
	  ar & m_vFloats;
	  ar & m_vmiType;
	  ar & m_vSeq;
	  ar & m_dIsotopeCorrection;
	  ar & m_iWindowCount;
	  ar & m_iConvolutionWidth;
	  ar & m_maxEnd;
	  ar & m_bDebug;
	  // omitted	double *m_pFactorial;
  }
  class vFloat : public std::vector<float> {
  public:
	    template<class Archive>
		void serialize(Archive &ar, const unsigned int version)
		{
			ar & *this;
		}
        float& at(int i)
        {
            return (*this)[i];
        }
        void init(int newSize)
        {
            resize(newSize);
            memset(&((*this)[0]), 0, newSize * sizeof(float));
        }
  };
#else
    class vFloat // ???? Replace this with stl vector
    {
    public:
        vFloat(void)
        {
            _size = _capacity = 0;
            _pf = NULL;
        }

        vFloat(const vFloat& that)
        {
            _size = _capacity = that._size;
            if (_size > 0)
            {
                _pf = (float *) malloc(_size * sizeof(float));
                memcpy(_pf, that._pf, _size * sizeof(float));
            }
            else
            {
                _pf = NULL;
            }
        }

        ~vFloat(void)
        {
            if (_pf != NULL)
                free(_pf);
        }

        void clear(void)
        {
            _size = 0;
        }

        int size(void)
        {
            return _size;
        }

        void ensure(int newSize)
        {
            if (newSize <= _capacity)
                return;

            if (_pf != NULL)
                free(_pf);
            _pf = (float *) malloc(newSize * sizeof(float));
            _capacity = newSize;
        }

        void init(int newSize)
        {
            ensure(newSize);
            _size = newSize;
            if (_pf == NULL)
            {
                assert(_size == 0);
                return;
            }
            memset(_pf, 0, _size * sizeof(float));
        }

        vFloat& operator=(const vFloat& that)
        {
            ensure(that._size);
            _size = that._size;
            if (_pf == NULL || that._pf == NULL)
            {
                assert(_size == 0);
                return *this;
            }
            memcpy(_pf, that._pf, that._size * sizeof(float));
            return *this;
        }

        float& operator[](int i)
        {
            return _pf[i];
        }

        float& at(int i)
        {
            assert(i >= 0 && i < _size);
            return _pf[i];
        }

    private:
        int _capacity;
        int _size;
        float *_pf;
    };
#endif // HAVE_MULTINODE_TANDEM

protected:
    friend class mscorefactory_c;
    mscore_c(void); // Should only be created by corresponding mscorefactory

public:
    virtual ~mscore_c(void);

    // Pluggable scoring API overrides
    virtual bool load_param(XmlParameter &_x);
    virtual void prescore(const size_t _i);
    virtual bool precondition(mspectrum &_s);
    virtual bool load_seq(const unsigned long _t, const long _c);
    virtual bool add_mi(mspectrum &_s);
    virtual double sfactor();
		virtual double hfactor(long _l); // hyper scoring factor from number of ions matched
    virtual unsigned long mconvert(double _m, const long _c);
    virtual void report_score(char *_buff, float _h);
    virtual bool clear();
		bool add_mi_hr(mspectrum &_s);
		virtual float score(const size_t _i); //MH: Checking other functions for scoring
		virtual float hconvert(float _h); // convert hyper score to histogram score

protected:
    bool postprocess_seq(const unsigned long _t, const long _c);
		bool postprocess_seq_hr(const unsigned long _t, const long _c);
    virtual double dot(unsigned long *_v);
		double dot_hr(unsigned long *_v);

    // ???? move this to mscore.h to be used by k-score, et al
    unsigned long imass(double _m)
    {
        return (unsigned long) (_m / m_dIsotopeCorrection + 0.5);
    }

    int neutralLoss(const int m, const float delta, const long z);

    void normalize_window(vFloat& processed,
                          float fWindowMaxI, float fMinCutoff,
                          vmiIter itWindowMI, vmiIter itWindowEnd);

    vFloat m_miUsed;
    vector<vFloat> m_vFloats;

		//MH: These are for accurate mass
		vectorvmiType m_vmiType;

    vmiType m_vSeq;

    double m_dIsotopeCorrection;
    int m_iWindowCount;
    int m_iConvolutionWidth;
    int m_maxEnd;
    bool m_bDebug;

		double *m_pFactorial;

    static bool byMassAsc(const MIType &_l, const MIType &_r)
    {
        return _l.m_lM < _r.m_lM;
    }
};

/*
 * mscorefactory_c implements a factory for creating mscore_c instances.
 */
class mscorefactory_c : public mpluginfactory
{
public:
    mscorefactory_c();

    virtual mplugin* create_plugin();
};
#endif // MSCORE_C_H
