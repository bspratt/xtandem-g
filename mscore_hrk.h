/*
 * Portions are Copyright (c) 2003-2006 Fred Hutchinson Cancer Research Center
 * Additional code Copyright (c) 2010-2011 Institute for Systems Biology
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
#include "mscore_k.h"

class mscore_hrk : public mscore
{
#ifdef HAVE_MULTINODE_TANDEM
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
      ar & boost::serialization::base_object<mscore>(*this);
      ar & m_maxEnd;
      ar & m_miUsed;
      ar & m_vmiType;
      ar & m_vmiUsed;
      ar & m_dIsotopeCorrection;
      
  }
#endif
protected:
    friend class mscorefactory_hrk;

    mscore_hrk(void);    // Should only be created through mscorefactory_tandem

public:
    virtual ~mscore_hrk(void);

public:
		virtual bool load_param(XmlParameter &_x); // allows score object to issue warnings,
                                               // or set variables based on xml
    virtual bool precondition(mspectrum &_s); // called before spectrum conditioning
    virtual void prescore(const size_t _i); // called before scoring
    
    virtual bool add_mi(mspectrum &_s);
		bool add_mi_hr(mspectrum &_s);

		virtual double sfactor(); // factor applied to final convolution score
    virtual unsigned long mconvert(double _m, const long _c); // convert mass to integer ion m/z for mi vector
    virtual void report_score(char* _buff, float _h); // format hyper score for output

    virtual bool clear();

protected:
    virtual double dot(unsigned long *_v); // this is where the real scoring happens
		double dot_hr(unsigned long *_v);
		virtual bool add_A(const unsigned long _t,const long _c);
		virtual bool add_B(const unsigned long _t,const long _c);
		virtual bool add_C(const unsigned long _t,const long _c);
		virtual bool add_Y(const unsigned long _t,const long _c);
		virtual bool add_X(const unsigned long _t,const long _c);
		virtual bool add_Z(const unsigned long _t,const long _c);

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
 * mscorefactory_hrk implements a factory for creating mscore_hrk instances.
 */
class mscorefactory_hrk : public mpluginfactory
{
public:
    mscorefactory_hrk();

    virtual mplugin* create_plugin();
};
