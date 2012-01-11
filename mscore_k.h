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
#include "mscore.h"

#ifndef MSCORE_K_H
#define MSCORE_K_H

class miLookup
{
#ifdef HAVE_MULTINODE_TANDEM // support for Hadoop and/or MPI?
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
		int initial_len = m_end - m_start;
        ar & m_start;
        ar & m_end;
		int len = m_end - m_start;
		if (len > initial_len) { // we're loading
			init(m_start, m_end);
		}
		while (len--) {
			ar & m_pfI[len];
		}
    }
#endif // HAVE_MULTINODE_TANDEM
public:
    miLookup(void) {
        init();
    }

    miLookup(const miLookup& rhs) {
        init();
        (*this) = rhs;
    }

    void init() {
        m_start = 0;
        m_end = 0;
        m_pfI = NULL;
    }

    virtual ~miLookup(void)    {
        if (m_pfI != NULL)
            free(m_pfI);
    }

    void init(int start, int end) {
        m_start = start;
        m_end = end;
        if (m_pfI != NULL)
            free(m_pfI);
        m_pfI = (float*) calloc(end - start, sizeof(float));
    }

    void clear(void) {
        memset(m_pfI, 0, (m_end - m_start) * sizeof(float));
    }

    int m_start; // start dalton value of pfI[0]
    int m_end; // end dalton value after pfI (i.e. end - start == size)
    float* m_pfI; // intensity lookup for masses rounded to dalton int values

    float& operator[](int i) {
        return m_pfI[i - m_start];
    }

    float get(int i) {
        if(i < m_start || i >= m_end)
            return 0.0;
        return m_pfI[i - m_start];
    }

/*
 * a simple copy operation, using the = operator
 */
    miLookup& operator=(const miLookup &rhs)    {
        m_start = rhs.m_start;
        m_end = rhs.m_end;

        if (m_pfI != NULL)
            free(m_pfI);
        size_t size = (m_end - m_start) * sizeof(float);
        m_pfI = (float*) malloc(size);
        memcpy(m_pfI, rhs.m_pfI, size);

        return *this;
    }
};

class mscore_k : public mscore
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
    friend class mscorefactory_k;

    mscore_k(void);    // Should only be created through mscorefactory_tandem

public:
    virtual ~mscore_k(void);

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
 * mscorefactory_k implements a factory for creating mscore_k instances.
 */
class mscorefactory_k : public mpluginfactory
{
public:
    mscorefactory_k();

    virtual mplugin* create_plugin();
};

#endif // MSCORE_K_H
