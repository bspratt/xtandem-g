/*
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

#include "stdafx.h"
#include <float.h>
#include "msequence.h"
#include "mspectrum.h"
#include "msequtilities.h"
#include "xmlparameter.h"
#include "mscore_c.h"

// Factory instance, registers itself with the mscoremanager.
static mscorefactory_c factory;
	
mscorefactory_c::mscorefactory_c()
{
    mscoremanager::register_factory("c-score", this);
}

mplugin* mscorefactory_c::create_plugin()
{
    return new mscore_c();
}

mscore_c::mscore_c(void)
{
    m_maxEnd = 0;
    m_dScale = 1.0; //0.05;
    m_dIsotopeCorrection = 1.0;
    m_iWindowCount = 10;
    m_iConvolutionWidth = 75;

	//MH: Set up factorial
	m_pFactorial = new double[64];
	double dFac = 1.0;
	m_pFactorial[0] = 1.0;
	long a = 1;
	while(a < 64)	{
		dFac *= (double)a;
		m_pFactorial[a] = dFac;
		a++;
	}
}

mscore_c::~mscore_c(void)
{
}

void mscore_c::prescore(const size_t _i)
{
    mscore::prescore(_i);
    m_miUsed.init(m_maxEnd);
}

bool mscore_c::clear()
{
    // ???? Need to call superclass?
    // mscore::clear();
    m_vFloats.clear();
		m_vmiType.clear(); //MH: for high-res
    return true;
}

/*
 * allows score object to issue warnings or set variable based on xml.
 */
bool mscore_c::load_param(XmlParameter &_x)
{
    if (!mscore::load_param(_x))
    return false;

    if (m_pSeqUtilFrag == &m_seqUtil)
        m_dIsotopeCorrection = 1.0005; /* monoisotopic */
    else
        m_dIsotopeCorrection = 1.0011; /* average */

    string strValue;
    string strKey = "c-score, histogram scale";
    if(_x.get(strKey,strValue))	{
        m_dScale = atof(strValue.c_str());
    }

    strKey = "c-score, debug";
    m_bDebug = false;
    if(_x.get(strKey,strValue))	{
        m_bDebug = 0 == strValue.compare("yes");
        cerr << "Debug is set to " << m_bDebug << endl;
    }

    return true;
}

/*
 * called before spectrum conditioning to allow the score object to
 * modify the spectrum in ways specific to the scoring algorithm.
 * default implementation does nothing.
 */
bool mscore_c::precondition(mspectrum &_s)
{
    if (_s.m_vMI.size() == 0)
        return false;

    if (!mscore::precondition(_s))
        return false;

    return true;
}


/*
 *
 */
void mscore_c::normalize_window(vFloat& processed, float fWindowMaxI, float fMinCutoff,
                                vmiIter itWindowMI,
                                vmiIter itWindowEnd)
{
    fWindowMaxI = sqrt(fWindowMaxI);

    //if (fWindowMaxI <= fMinCutoff)
    //    return;

    float fScaleI = 50.f / fWindowMaxI;
    for ( ; itWindowMI != itWindowEnd; ++itWindowMI)
    {
        if (itWindowMI->m_fI > fMinCutoff)
        {
            int iM = imass(itWindowMI->m_fM);
            float fI = sqrt(itWindowMI->m_fI) * fScaleI;
            processed[iM] = max(processed[iM], fI);
        }
    }
}

/*
 * add_mi does the work necessary to set up an mspectrum object for modeling. 
 *   - an entry in the m_State object is made for the parent ion M+H
 * once an mspectrum has been added, the original mspectrum is no longer
 * needed for modeling, as all of the work associated with a spectrum
 * is only done once, prior to modeling sequences.
 */
bool mscore_c::add_mi(mspectrum &_s)
{
	//MH - intercept function for high res data
	if(m_lErrorType & mscore::T_FRAGMENT_PPM){
		return add_mi_hr(_s);
	} else if ( (m_lErrorType & mscore::T_FRAGMENT_DALTONS) && m_fErr<1.0){

		//Must request smaller than 1 dalton bin to use high-res
		return add_mi_hr(_s);
	}

    const float fThresholdI = 0.05f; // 5% max intensity cut-off

    if (!mscore::add_mi(_s))
        return false;

    vFloat processed;

    if (_s.m_vMI.size() == 0)
    {
        m_vFloats.push_back(processed);
        return true;
    }

    vmiIter itMI = _s.m_vMI.begin();
    vmiIter itEnd = _s.m_vMI.end();

    int startMass =  imass(itMI->m_fM);
    int endMass = imass(itEnd[-1].m_fM);
    int size = endMass + 2 * m_iConvolutionWidth;

    if (size > m_maxEnd)
        m_maxEnd = size;

    float fMaxI = 0.f;
    for (itMI = _s.m_vMI.begin(); itMI != itEnd; ++itMI)
    {
        float fI = itMI->m_fI;
        if (fI > fMaxI)
            fMaxI = fI;
    }
    float fMinCutoff = fThresholdI * sqrt(fMaxI);

    processed.init(size);

    float fWindowSize = endMass * 1.f / m_iWindowCount;

    int iWindowIdx;
    int iWindowEnd;
    float fWindowMaxI = 0.f;
    vmiIter itWindowMI = _s.m_vMI.begin();

    for (int i = 0; i < m_iWindowCount; ++i)
    {
        iWindowEnd = (int) ((i + 1) * fWindowSize);
        if (iWindowEnd > startMass)
        {
            iWindowIdx = i;
            break;
        }
    }

    for (itMI = _s.m_vMI.begin(); itMI != itEnd; ++itMI)
    {
        int iM = imass(itMI->m_fM);
        if (iM < iWindowEnd)
        {
            if (itMI->m_fI > fWindowMaxI)
                fWindowMaxI = itMI->m_fI;
        }
        else
        {
            // normalize current window
            normalize_window(processed, fWindowMaxI, fMinCutoff, itWindowMI, itMI);

            // update params for next window
            ++iWindowIdx;
            iWindowEnd = (int) ((iWindowIdx + 1) * fWindowSize);
            fWindowMaxI = 0.f;
            itWindowMI = itMI;
        }
    }

    // normalize last window
    normalize_window(processed, fWindowMaxI, fMinCutoff, itWindowMI, itEnd);

    double frac = 1.0 / (2.0 * m_iConvolutionWidth + 1.0);

    vector<float> tempRangeLookup;
    tempRangeLookup.resize(size, 0.f);

    for (int i = max(0, startMass - m_iConvolutionWidth); i < endMass + m_iConvolutionWidth; ++i)
    {
        double sum = 0.0;
        for (int ii = max(0, i - m_iConvolutionWidth); ii <= i + m_iConvolutionWidth ; ++ii)
            sum += processed[ii];
        tempRangeLookup[i] = (float) (sum * frac);
    }

    for (int i = 0; i < size; ++i)
    {
        processed[i] -= tempRangeLookup[i];
#ifdef DEBUGABLE
        if (m_bDebug)
            cerr << "Processed\t" << i << "\t" << processed[i] << endl;
#endif
    }

    m_vFloats.push_back(processed);

    return true;
}

/*
 * mconvert converts from mass and charge to integer ion value
 * for mi vector.
 */
unsigned long mscore_c::mconvert(double _m, const long _c)
{
    double dMass = _m / _c + m_pSeqUtilFrag->m_dProton;
    return imass(dMass);
}

/*
 * sfactor returns a factor applied to the final convolution score.
 */
double mscore_c::sfactor()
{
    // ???? insert length dependence here?
    return .02;
}

/*
 * report_score formats a hyper score for output.
 */
void mscore_c::report_score(char* buffer, float hyper)
{
    //sprintf(buffer, "%d",(int) (hyper + 0.5f));
	sprintf(buffer,"%.1f",hconvert(hyper));
}

__inline__ int mscore_c::neutralLoss(const int m, const float delta, const long z)
{
    return (int) (m - delta / z + 0.5f);
}

bool mscore_c::postprocess_seq(const unsigned long _t, const long _c)
{
	//MH - intercept function for high res data
	if(m_lErrorType & mscore::T_FRAGMENT_PPM){
		return postprocess_seq_hr(_t,_c);
	} else if ( (m_lErrorType & mscore::T_FRAGMENT_DALTONS) && m_fErr<1.0){

		//Must request smaller than 1 dalton bin to use high-res
		return postprocess_seq_hr(_t,_c);
	}

    m_vSeq.clear();
    MIType uSeq;

    for (int a = 0; m_plSeq[a] != 0; ++a) {
        int i = (int) m_plSeq[a];

        uSeq.m_fI = 50.f;
        uSeq.m_lM = i;
        m_vSeq.push_back(uSeq);

        uSeq.m_fI = 25.f;
        uSeq.m_lM = i - 1;
        m_vSeq.push_back(uSeq);
        uSeq.m_lM = i + 1;
        m_vSeq.push_back(uSeq);

        // ???? parameter that will allow the user to exclude even for _c==2?
        // if charge is >= 3, don't clutter things up with neutral loss peaks
        if (_c < 3 && (T_B | T_Y) & _t)
        {
            uSeq.m_fI = 10.f;
            uSeq.m_lM = neutralLoss(i, 17.f, _c);
            m_vSeq.push_back(uSeq);
            uSeq.m_lM = neutralLoss(i, 18.f, _c);
            m_vSeq.push_back(uSeq);
            if (T_B & _t)
            {
                uSeq.m_lM = neutralLoss(i, 28.f, _c);
                m_vSeq.push_back(uSeq);
            }
        }
    }

    return true;
}

bool mscore_c::load_seq(const unsigned long _t,const long _c)
{
    if (!mscore::load_seq(_t, _c))
        return false;
    return postprocess_seq(_t, _c);
}

/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses).
 */
double mscore_c::dot(unsigned long *_v)
{
	//MH - intercept function for high res data
	if(m_lErrorType & mscore::T_FRAGMENT_PPM){
		return dot_hr(_v);
	} else if ( (m_lErrorType & mscore::T_FRAGMENT_DALTONS) && m_fErr<1.0){

		//Must request smaller than 1 dalton bin to use high-res
		return dot_hr(_v);
	}

    unsigned long lCount = 0;

    vFloat& processed = m_vFloats[m_lId];

    int size = processed.size();

    vmiType::iterator itFrag = m_vSeq.begin();
    vmiType::const_iterator itFragEnd = m_vSeq.end();

    double dScore = 0.0;

    for (; itFrag != itFragEnd; ++itFrag)
    {
        int iIon = (int) itFrag->m_lM;
        float fIntensity = itFrag->m_fI;

        if (iIon >= size)
            continue; // break;

        float fTmp = processed[iIon] * fIntensity;
        float previous = m_miUsed[iIon];

        if ((fTmp != 0.f && previous == 0.f) || previous < fTmp) {
            dScore += fTmp - previous;
            m_miUsed[iIon] = fTmp;
            if (fIntensity >= 50)
                ++lCount; // hack; just count primary peak matches
        }
    }
    *_v = lCount;
    return (dScore);
}

static void dumpFrag(vmiType& vSeq)
{
    vector<MIType>::iterator itFrag = vSeq.begin();
    vector<MIType>::const_iterator itFragEnd = vSeq.end();

    for (; itFrag != itFragEnd; ++itFrag)
        printf("%d %f\n", itFrag->m_lM, itFrag->m_fI);
}

bool mscore_c::add_mi_hr(mspectrum &_s)
{
    const float fThresholdI = 0.05f; // 5% max intensity cut-off

		vmiType vType;

    if (!mscore::add_mi(_s)) return false;

    vFloat processed;

    if (_s.m_vMI.size() == 0) {
			m_vmiType.push_back(vType);
      return true;
    }

		int i;
		int size=(int)_s.m_vMI.size();
		float startMass = _s.m_vMI[0].m_fM;
		float endMass = _s.m_vMI[size-1].m_fM;

		//See if we have a new maximum spectrum size
		//Note that this is always an overestimate assuming at least one peak
		//will not pass the above threshold
		if (size > m_maxEnd) m_maxEnd = size;

		//Find the spectrum max
    float fMaxI = 0.f;
    for (i = 0; i < size; i++) {
			float fI = _s.m_vMI[i].m_fI;
			if (fI > fMaxI) fMaxI = fI;
    }
    float fMinCutoff = fThresholdI * sqrt(fMaxI);

    processed.init(size);

    float fWindowSize = endMass / m_iWindowCount;

    int iWindowIdx=1;
    float fWindowEnd;
    float fWindowMaxI = 0.f;

		//Normalize across each window
		int j=0;
		int iFirst=0;
		fWindowEnd=fWindowSize;
		for(j=0;j<size;j++){
			float fM = _s.m_vMI[j].m_fM;
      if (fM < fWindowEnd) {
				if (_s.m_vMI[j].m_fI > fWindowMaxI) fWindowMaxI = _s.m_vMI[j].m_fI;
			} else {
				fWindowMaxI=sqrt(fWindowMaxI);
				for(i=iFirst;i<j;i++){
					if (_s.m_vMI[i].m_fI > fMinCutoff) {
						processed[i]= sqrt(_s.m_vMI[i].m_fI)*50.0f/fWindowMaxI;
					}
				}
				iFirst=j;
				iWindowIdx++;
				fWindowEnd = iWindowIdx* fWindowSize;
        fWindowMaxI = _s.m_vMI[j].m_fI;
			}
		}

    // normalize last window
    fWindowMaxI=sqrt(fWindowMaxI);
		for(i=iFirst;i<j;i++) {
			if (_s.m_vMI[i].m_fI > fMinCutoff){
				processed[i]= sqrt(_s.m_vMI[i].m_fI)*50.0f/fWindowMaxI;
			}
		}

    double frac = (2.0 * m_iConvolutionWidth + 1.0);

    vector<float> tempRangeLookup;
    tempRangeLookup.resize(size, 0.f);

		//Use that funky thing
		for(i=0;i<size;i++){
			double sum=0.0;
			for(j=0; j<size && _s.m_vMI[j].m_fM <= _s.m_vMI[i].m_fM + m_iConvolutionWidth;j++){
				if(_s.m_vMI[j].m_fM < _s.m_vMI[i].m_fM - m_iConvolutionWidth) continue;
				sum+=processed[j];
			}
			tempRangeLookup[i] = (float) (sum/frac);
		}
		
		//Make array of normalized datapoints.
		//Multiply by 1000 to suppress loss of precision due to using integers
		//Using integers prevents creating yet another data type.
		MIType uType;
		for(i=0;i<size;i++) {
			if(processed[i]>tempRangeLookup[i]){
				uType.m_lM = (int)((processed[i]-tempRangeLookup[i])*100.0f);
				uType.m_fI = _s.m_vMI[i].m_fM;
				vType.push_back(uType);
			}
		}

    m_vmiType.push_back(vType);

    return true;
}

double mscore_c::dot_hr(unsigned long *_v)
{
  unsigned long lCount = 0;
	int iTmp;
	double dScore = 0.0;

	float fPpm;
	float fPpmUser;
	float fPpmUser2x;

  vmiType& processed = m_vmiType[m_lId];

  int size = processed.size();
	if (size<2){
		_v=0;
		return dScore;
	}

	/*
	for(int j=0;j<size;j++){
		cout << processed[j].m_fI << " " << processed[j].m_lM << endl;
	}
	for(int j=0;j<m_vSeq.size();j++){
		cout << m_vSeq[j].m_fI << " " << m_vSeq[j].m_lM << endl;
	}
	*/

	//Treat tolerances different depending on PPM or Daltons
	if(m_lErrorType & T_FRAGMENT_PPM) fPpmUser = (float)(m_fErr*1e6/200.0); //convert back to true ppm parameter
	else fPpmUser = m_fErr;

	//Check out twice as far for half the score - this is equivalent to checking the next bin
	fPpmUser2x = (float)(fPpmUser*2.0);

	int i=1;
	int iBin;
	int iSeqSize=m_vSeq.size();
	
	//perform a single pass through each array.
	//check every point in m_pfSeq, but don't revisit positions in m_vmiType
	for (int a = 0; a < iSeqSize; a++) {

		float fIon = m_vSeq[a].m_fI;
		int iIntensity = (int)m_vSeq[a].m_lM;

		while(fIon > processed[i].m_fI){
			i++;
			if(i==size){
				i--;
				break;
			}
		}
		if(fIon > processed[i].m_fI+1.0) break;

		//Use different calculation based on unit type
		if( (processed[i].m_fI-fIon) < (fIon-processed[i-1].m_fI) ){
			if(m_lErrorType & T_FRAGMENT_PPM)	fPpm = (float)(-(fIon-processed[i].m_fI)/processed[i].m_fI*1e6);
			else fPpm = processed[i].m_fI - fIon;
			iBin=i;
		} else {
			if(m_lErrorType & T_FRAGMENT_PPM) fPpm = (float)((fIon-processed[i-1].m_fI)/processed[i-1].m_fI*1e6);
			else fPpm = fIon - processed[i-1].m_fI;
			iBin=i-1;
		}
		//for some cases, it is possible to still have a negative ppm
		if(fPpm<0.0f) fPpm=-fPpm;

		//cout << fIon << "\t" << processed[iBin].m_fI << "\t" << fPpm << endl;

		//Check within tolerance
		if(fPpm<fPpmUser){
			//cout << "wtf: " << fIon << "\t" << processed[iBin].m_fI << "\t" << fPpm << endl;
			iTmp = processed[iBin].m_lM*iIntensity;
      if ((int)m_miUsed[iBin] < iTmp) {
        dScore += (double)(iTmp - (int)m_miUsed[iBin]);
        m_miUsed[iBin] = (float)iTmp;
        if(iIntensity==50) lCount++;
      }
			continue;
    }

		//Check at twice the tolerance
		if(fPpm<fPpmUser2x){
			//cout << "wtfx2: " << fIon << "\t" << processed[iBin].m_fI << "\t" << fPpm << endl;
			iTmp = processed[iBin].m_lM*iIntensity;
			if(iIntensity==50){
				iTmp /= 2;
				if ((int)m_miUsed[iBin] < iTmp) {
					dScore += (double)(iTmp - (int)m_miUsed[iBin]);
					m_miUsed[iBin] = (float)iTmp;
				}
			}
    }

	}

  *_v = lCount;
	//cout << lCount << "\t" << dScore/100.0 << endl;
  return (dScore/100.0);
}

bool mscore_c::postprocess_seq_hr(const unsigned long _t, const long _c)
{
    m_vSeq.clear();
    MIType uSeq;

    for (int a = 0; m_pfSeq[a] != 0; ++a) {
 
        uSeq.m_fI = m_pfSeq[a];
        uSeq.m_lM = 50;
        m_vSeq.push_back(uSeq);

        // ???? parameter that will allow the user to exclude even for _c==2?
        // if charge is >= 3, don't clutter things up with neutral loss peaks
        if (_c < 3 && (T_B | T_Y) & _t)
        {
            uSeq.m_lM = 10;
            uSeq.m_fI = m_pfSeq[a]-17.026549105f;
            m_vSeq.push_back(uSeq);
            uSeq.m_fI = m_pfSeq[a]-18.0105647f;
            m_vSeq.push_back(uSeq);
            if (T_B & _t)
            {
                uSeq.m_fI = m_pfSeq[a]-27.99491463f;
                m_vSeq.push_back(uSeq);
            }
        }
    }

    return true;
}

/*
 * the score method is called externally to score a loaded peptide against one of the
 * loaded mass spectra. the mass spectrum is refered to by its index number in the
 * m_vSpec mspectrum vector. the sequence has already been loaded via set_seq and/or
 * add_seq.
 */
float mscore_c::score(const size_t _i)
{
	m_fScore = -1.0;
	m_fHyper = -1.0;
	double dFactor = 1.0;
/*
 * return -1000.0 if there is no sequence available
 */
	if(m_pSeq == NULL)
		return -1000.0;
/*
 * initialize values for the protein modeling session
 */
	prescore(_i);

	double dScore = (float)0.0;
	double dValue = (float)0.0;

	unsigned long lType = T_Y;
	unsigned long lValue = 0;
	unsigned long lValueTotal = 0;
	unsigned long lS = S_Y;

	//MH: I'm not sure I like this bit.
	long lChargeLimit = (long)m_vSpec[m_lId].m_fZ;
	if(lChargeLimit == 1)	{
		lChargeLimit = 2;
	}
	if((m_lType & T_C) || (m_lType & T_Z))	{
		if(lChargeLimit > 2)	{
			lChargeLimit--;
		}
	}

/*
 * iterate through all of the possible values of the mscore_type enum
 * comparing them against m_lType
 */
	while(lType < m_lType+1)	{
		lValueTotal = 0;
		dValue = 0.0;
		if(lType & m_lType)	{
			long a = 1;
			while(a < lChargeLimit)	{
/*
 * load the sequence arrays for each possible charge states for the selected spectrum
 */
				load_seq(lType,a);
				lValue = 0;
/*
 * perform a dot product on each charge state
 */
				dValue += dot(&lValue);

				//MH: I'm not sure why there is this special case for y-ions in the first charge series of
				//a plus two precursor
				/*
				if(a == 1 && T_Y & lType && (long)m_vSpec[m_lId].m_fZ == 2)	{
					unsigned long lTemp = 0;
					add_Y(0,2);
					dValue += dot(&lTemp);
					lValue += lTemp;
				}
				*/
				lValueTotal += lValue;
				a++;
			}
			dScore += dValue;
		}
		else m_pfScore[lS] = (float) (dValue*sfactor());
		m_plCount[lS] = lValueTotal;
/*
 * move on to next value in the mstate_type enum
 */
		lS++;
		lType *= 2;
	}

	dScore *= sfactor();
	m_fScore = (float)dScore;
	
	//MH: attempts to replicate orignal tandem scoring escalator
	dFactor = dScore;
	lType = T_Y;
	lS = S_Y;
	long lNScale=0;
	long lCScale=0;
	while(lType < m_lType+1)	{
		//if(lType & m_lType)	{
		//	dFactor*= hfactor(m_plCount[lS]);
		//}
		if(lType & m_lType && (lType & T_X || lType & T_Y || lType & T_Z) )	lCScale+=m_plCount[lS];
		if(lType & m_lType && (lType & T_A || lType & T_B || lType & T_C) )	lNScale+=m_plCount[lS];
		
		lS++;
		lType*=2;
	}
	
	dFactor*=hfactor(lNScale);
	dFactor*=hfactor(lCScale);

	if(dFactor > FLT_MAX) m_fHyper = (float)FLT_MAX;
	else m_fHyper = (float)dFactor;
	if(m_fHyper<1.0f) m_fHyper=1.0f;

	//if(dFactor > FLT_MAX) m_fHyper = (float)log10(FLT_MAX);
	//else m_fHyper = (float)log10(dFactor);
	//if(m_fHyper<0.0f) m_fHyper=0.0f;

/*
 * returning 1.0 for a zero score makes the logic in mprocess easier. see mprocess:create_score
 * to see why.
 */
	if(dScore == 0.0)	{
		dScore = 1.0;
		m_fHyper = 1.0f;
	}
	
	//cout << dScore << "\t" << m_fHyper << "\t" << dFactor << endl;
	return (float) dScore;
}

double mscore_c::hfactor(long _l) {
	return m_pFactorial[_l];
}

float mscore_c::hconvert(float _f) {
	if(_f <= 0.0)
		return 0.0;
	return (float)(m_dScale*log10(_f));
	//return (float)(m_dScale*_f);
}