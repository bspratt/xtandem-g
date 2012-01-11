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

#include "stdafx.h"
#include "msequence.h"
#include "mspectrum.h"
#include "msequtilities.h"
#include "xmlparameter.h"
#include "mscore_hrk.h"

// Factory instance, registers itself with the mscoremanager.
static mscorefactory_hrk factory;
	
mscorefactory_hrk::mscorefactory_hrk()
{
    mscoremanager::register_factory("hrk-score", this);
}

mplugin* mscorefactory_hrk::create_plugin()
{
    return new mscore_hrk();
}

mscore_hrk::mscore_hrk(void)
{
    m_dScale = 0.05;
    m_maxEnd = 0;
    m_dIsotopeCorrection = 1.0;
}

mscore_hrk::~mscore_hrk(void)
{
}

bool mscore_hrk::clear()
{
    m_vmiType.clear();
    return true;
}

/*
 * allows score object to issue warnings, or set variable based on xml.
 */
bool mscore_hrk::load_param(XmlParameter &_x)
{
    if (!mscore::load_param(_x))
    return false;

    if (m_pSeqUtilFrag == &m_seqUtil)
        m_dIsotopeCorrection = 1.0005; /* monoisotopic */
    else
        m_dIsotopeCorrection = 1.0011; /* average */

    string strValue;
    string strKey = "k-score, histogram scale";
    if(_x.get(strKey,strValue))	{
        m_dScale = atof(strValue.c_str());
    }

    return true;
}

/*
 * called before spectrum conditioning to allow the score object to
 * modify the spectrum in ways specific to the scoring algorithm.
 * default implementation does nothing.
 */
bool mscore_hrk::precondition(mspectrum &_s)
{
    if (_s.m_vMI.size() == 0)
        return false;

    if (!mscore::precondition(_s))
        return false;

    return true;
}

/*
 * called before scoring inside the score() function to allow any
 * necessary resetting of member variables.
 */
void mscore_hrk::prescore(const size_t _i)
{
    mscore::prescore(_i);

    // Initialize of clear the used intensity look-up.
    if (m_miUsed.m_pfI == NULL)
        m_miUsed.init(0, m_maxEnd);
    else
        m_miUsed.clear();
}

/*
 * mconvert converts from mass and charge to integer ion value
 * for mi vector.
 */
unsigned long mscore_hrk::mconvert(double _m, const long _c)
{
    const double fZ = (double)_c;

    double dMass = (fZ*m_pSeqUtilFrag->m_dProton + _m)/fZ;
    return imass(dMass);
}

/*
 * sfactor returns a factor applied to the final convolution score.
 */
double mscore_hrk::sfactor()
{
/*
 * Multiply by log(length) to remove length dependence on dot product score
 * Divide by 3.0 to scale score to 1000
 */
    double dFactor = log((double)m_lSeqLength)*1.0 /
        (3.0*sqrt((double)m_lSeqLength)); /* change iLenPeptide to tot # of fragment ions? */
    dFactor *= 1000.0;
    return dFactor;
}

/*
 * report_score formats a hyper score for output.
 */
void mscore_hrk::report_score(char* buffer, float hyperscore)
{
    sprintf(buffer, "%d",(int) (hyperscore + 0.5));
}

/*
 * add_mi does the work necessary to set up an mspectrum object for modeling. 
 *   - an entry in the m_State object is made for the parent ion M+H
 * once an mspectrum has been added, the original mspectrum is no longer
 * needed for modeling, as all of the work associated with a spectrum
 * is only done once, prior to modeling sequences.
 */
bool mscore_hrk::add_mi(mspectrum &_s)
{
    if (!mscore::add_mi(_s))
        return false;

		//MH - intercept function for high res data
		if(m_lErrorType & mscore::T_FRAGMENT_PPM){
			return add_mi_hr(_s);
		} else if ( (m_lErrorType & mscore::T_FRAGMENT_DALTONS) && m_fErr<1.0){

			//Must request smaller than 1 dalton bin to use high-res
			return add_mi_hr(_s);
		}

    vmiType vType;
    if (_s.m_vMI.size() == 0)
    {
        // Return value appears to be ignored, so just add empty type.
        m_vmiType.push_back(vType);
        return true;
    }

    int iWindowCount = 10;
    float fMaxI = 0;
    float fTotI = 0;

    miLookup tempLookup;

    vector<mi>::iterator itMI = _s.m_vMI.begin();
    vector<mi>::iterator itEnd = _s.m_vMI.end();
    int startMass = imass(itMI->m_fM);
    int endMass = imass(itEnd[-1].m_fM);

    // Screen peeks on upper end.
    int endMassMax = (int)(((_s.m_dMH + (_s.m_fZ - 1) * m_seqUtil.m_dProton) / _s.m_fZ) * 2.0 + 0.5) + iWindowCount;
    while (itMI != itEnd && endMass >= endMassMax) {
        itEnd--;
        endMass = imass(itEnd[-1].m_fM);
    }

    if (itMI == itEnd)  // No peaks left.
    {
        // Return value appears to be ignored, so just add empty type.
        m_vmiType.push_back(vType);
        return true;
    }

    tempLookup.init(max(0, startMass - 50), endMass + 50);

    if (tempLookup.m_end > m_maxEnd)
        m_maxEnd = tempLookup.m_end;

    int peakCount = 0;
    while (itMI != itEnd) {

        int iM = imass(itMI->m_fM);
        float fI = sqrt(itMI->m_fI);

        fTotI += fI;
        if (fMaxI < fI)
            fMaxI = fI;

        if (tempLookup[iM] < fI)
        {
            if (tempLookup[iM] == 0.0)
                peakCount++;
            tempLookup[iM] = fI;
        }

        itMI++;
    }

    float fMinCutoff = (float) (0.05 * fMaxI);

    int range = (int) min(endMassMax, iWindowCount + endMass) - startMass;
    if (range > 3000)
        iWindowCount=10;
    else if (range > 2500)
        iWindowCount=9;
    else if (range > 2000)
        iWindowCount=8;
    else if (range > 1500)
        iWindowCount=7;
    else if (range > 1000)
        iWindowCount=6;
    else
        iWindowCount=5;

    int iWindowSize = range / iWindowCount;

    /*
     * Process input spectrum for dot product - split windows
     */
    for (int i = 0; i < iWindowCount; i++) {

        float fMaxWindowI = 0.0;
        int iStart = startMass + i*iWindowSize;

        /*
         * Get maximum intensity within window
         */
        for (int ii = iStart; ii < iStart + iWindowSize; ii++) {
            float fI = tempLookup[ii];
            if (fI > fMaxWindowI)
                fMaxWindowI = fI;
        }

        if (fMaxWindowI > 0.0 && fMaxWindowI > fMinCutoff) {
            double dFactor= 1.0 / fMaxWindowI;

            /*
             * Normalize within window
             */
            for (int ii = iStart; ii < iStart + iWindowSize; ii++) {
                double dI = tempLookup.get(ii);
                tempLookup[ii] = (float) (dI * fMaxI * dFactor);
            }
        }
    }

    /*
     * Reduce intensity and make unit vector by dividing
     * every point by sqrt(sum(x^2))
     */
    double dSpectrumArea = 0.0;
    for (int i = startMass; i <= endMass; i++) {
        double d = tempLookup.get(i);
        if (d > 0.0)
            dSpectrumArea += d * d;
    }

    dSpectrumArea = sqrt(dSpectrumArea);

    for (int i = startMass; i <= endMass; i++) {
        float f = tempLookup.get(i);
        if (f > 0.0)
            tempLookup[i] = (float) (f / dSpectrumArea);
    }

    /*
     * Perform mix-range modification to input spectrum
     */
    miLookup tempRangeLookup;
    tempRangeLookup.init(tempLookup.m_start, tempLookup.m_end);
    for (int i = tempLookup.m_start; i < tempLookup.m_end; i++) {
        double sum = 0.0;
        for (int ii = i - 50 ; ii <= i + 50 ; ii++)
            sum += tempLookup.get(ii);
        tempRangeLookup[i] = (float) (sum / 101.0);
    }

    MIType uType;
    for (int i = tempLookup.m_start; i < tempLookup.m_end; i++) {
        tempLookup[i] -= tempRangeLookup[i];
        if (tempLookup[i] > 0) {
            uType.m_lM = i;
            uType.m_fI = tempLookup[i];
            vType.push_back(uType);
        }
    }

    m_vmiType.push_back(vType);
    return true;
}

/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses).
 */
double mscore_hrk::dot(unsigned long *_v)
{

	//MH - intercept function for high res data
	if(m_lErrorType & mscore::T_FRAGMENT_PPM){
		return (dot_hr(_v));
	} else if ( (m_lErrorType & mscore::T_FRAGMENT_DALTONS) && m_fErr<1.0){

		//Must request smaller than 1 dalton bin to use high-res
		return (dot_hr(_v));
	}

    unsigned long lCount = 0;
    double dScore = 0.0;
    float fTmp;

    vector<MIType>::iterator itType = m_vmiType[m_lId].begin();
    vector<MIType>::const_iterator itEnd = m_vmiType[m_lId].end();
    vector<MIType>::const_iterator itBegin = itType;
    vector<MIType>::iterator itGreater;

		for (int a = 0; m_plSeq[a] != 0; a++) {

				int iIon = (int) m_plSeq[a];

				// Search for first peak in spectrum greater than or equal
				// that in current sequence.  Because there are usually a
				// lot more spectrum peaks than sequence peaks, jumping ahead
				// has significant performance benefits.

				const int step = 5;
				while (step < itEnd - itType && itType[step].m_lM < m_plSeq[a]) {
						itType += step;
				}
				while(itType != itEnd && itType->m_lM < m_plSeq[a]) {
						itType++;
				}

				itGreater = itType;
				if (itType != itEnd && iIon == itType->m_lM) {
						fTmp = itType->m_fI;
						if (m_miUsed.get(iIon) < fTmp) {
								dScore += fTmp - m_miUsed.get(iIon);
								m_miUsed[iIon] = fTmp;
								lCount++;
						}

						itGreater = itType + 1;
				}

				if (itType != itBegin && iIon - 1 == itType[-1].m_lM) {
						fTmp = ((float)0.5)*itType[-1].m_fI;
						if (m_miUsed.get(iIon-1) < fTmp) {
								dScore += fTmp - m_miUsed.get(iIon-1);
								m_miUsed[iIon-1] = fTmp;
						}
				}

				if (itGreater != itEnd && iIon + 1 == itGreater->m_lM) {
						fTmp = ((float)0.5)*itGreater->m_fI;
						if (m_miUsed.get(iIon+1) < fTmp) {
								dScore += fTmp - m_miUsed.get(iIon+1);
								m_miUsed[iIon+1] = fTmp;
						}
				}
		}

    *_v = lCount;    
    return (dScore);
}

/* MH
 * add_mi_hr is similar to add_mi for seting up an mspectrum object for modeling. 
 * It is designed for hi res data, and does not use integer binning. A paired
 * dot_hr function has been created to perform dot products from these models.
 *   - an entry in the m_State object is made for the parent ion M+H
 * once an mspectrum has been added, the original mspectrum is no longer
 * needed for modeling, as all of the work associated with a spectrum
 * is only done once, prior to modeling sequences.
 */
bool mscore_hrk::add_mi_hr(mspectrum &_s)
{

		vmiType vType;
    if (_s.m_vMI.size() == 0)
    {
        // Return value appears to be ignored, so just add empty type.
        m_vmiType.push_back(vType);
        return true;
    }

    int iWindowCount = 10;
    int iMaxI = 0;
    int iTotI = 0;

		MIType uType;
    vector<MIType> tempLookup;

    vector<mi>::iterator itMI = _s.m_vMI.begin();
    vector<mi>::iterator itEnd = _s.m_vMI.end();
    int startMass = imass(itMI->m_fM);
    int endMass = imass(itEnd[-1].m_fM);

    // Screen peeks on upper end.
    int endMassMax = (int)(((_s.m_dMH + (_s.m_fZ - 1) * m_seqUtil.m_dProton) / _s.m_fZ) * 2.0 + 0.5) + iWindowCount;
    while (itMI != itEnd && endMass >= endMassMax) {
        itEnd--;
        endMass = imass(itEnd[-1].m_fM);
    }

    if (itMI == itEnd)  // No peaks left.
    {
        // Return value appears to be ignored, so just add empty type.
        m_vmiType.push_back(vType);
        return true;
    }

    int peakCount = 0;
    while (itMI != itEnd) {

			  //note that these are reversed. Intensity is stored as an integer
			  //because it can be inaccurate. Mass is stored as a float for
			  //accuracy.
			  uType.m_lM = (int)(sqrt(itMI->m_fI)+0.5);
				uType.m_fI = itMI->m_fM;
				tempLookup.push_back(uType);

        if(iMaxI < (int)uType.m_lM) iMaxI = (int)uType.m_lM;

        itMI++;
    }

    float fMinCutoff = 0.05f * (float)iMaxI;

    int range = (int) min(endMassMax, iWindowCount + endMass) - startMass;
    if (range > 3000)
        iWindowCount=10;
    else if (range > 2500)
        iWindowCount=9;
    else if (range > 2000)
        iWindowCount=8;
    else if (range > 1500)
        iWindowCount=7;
    else if (range > 1000)
        iWindowCount=6;
    else
        iWindowCount=5;

    int iWindowSize = range / iWindowCount;

    /*
     * Process input spectrum for dot product - split windows
     */
		int i=0;
		int ii;
		int iStart=0;
		while(iStart<(int)tempLookup.size()){
			
			float fMaxWindowI = 0.0;
			float fWindowEnd = (float)startMass+i*iWindowSize;

			for(ii=iStart; ii<(int)tempLookup.size() && tempLookup[ii].m_fI<fWindowEnd; ii++){
				float fI = (float)tempLookup[ii].m_lM;
        if (fI > fMaxWindowI) fMaxWindowI = fI;
      }

			if (fMaxWindowI > 0.0 && fMaxWindowI > fMinCutoff) {

        //Normalize within window
        for (ii = iStart; ii<(int)tempLookup.size() && tempLookup[ii].m_fI<fWindowEnd; ii++) {
					double dI = (double)tempLookup[ii].m_lM;
					tempLookup[ii].m_lM = (int) (dI * iMaxI / fMaxWindowI);
        }
      }

			iStart=ii;
			i++;
		}

    /*
     * Reduce intensity and make unit vector by dividing
     * every point by sqrt(sum(x^2))
     */
    double dSpectrumArea = 0.0;
		for (i = 0; i < (int)tempLookup.size(); i++) {
			double d = (double) tempLookup[i].m_lM;
      if (d > 0.0) dSpectrumArea += d * d;
    }

    dSpectrumArea = sqrt(dSpectrumArea);

    for (i = 0; i < (int)tempLookup.size(); i++) {
			float f = (float)tempLookup[i].m_lM*1000.0f;
			if (f > 0.0) tempLookup[i].m_lM = (int) (f / dSpectrumArea);
    }

    /*
     * Perform mix-range modification to input spectrum
     */
    vector<int> tempRangeLookup;
    for (i = 0; i < (int)tempLookup.size(); i++) {
			double sum = 0;
			for (int ii = 0; tempLookup[ii].m_fI <= tempLookup[i].m_fI+50.0 ; ii++){
				if(tempLookup[ii].m_fI < tempLookup[i].m_fI-50.0) continue;
				sum += (double)tempLookup[ii].m_lM;
				if(ii==tempLookup.size()-1) break;
			}
			tempRangeLookup.push_back((int)(sum / 101.0));
    }

    for (i = 0; i < (int)tempLookup.size(); i++) {
			if( (int)tempLookup[i].m_lM > tempRangeLookup[i]) {
				tempLookup[i].m_lM -= tempRangeLookup[i];
				uType.m_lM = tempLookup[i].m_lM;
				uType.m_fI = tempLookup[i].m_fI;
				vType.push_back(uType);
      }
			if(i+1>m_maxEnd) m_maxEnd=i+1;
    }

    m_vmiType.push_back(vType);
    return true;
}

/* MH:
 * dot_hr is the same as dot, but modified for sparse arrays and high mass accuracy.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_pfSeq (float masses).
 */
double mscore_hrk::dot_hr(unsigned long *_v)
{

	unsigned long lCount = 0;
  double dScore = 0.0;
  int iTmp;

	float fPpm;
	float fPpmUser;
	float fPpmUser2x;

	//Treat tolerances different depending on PPM or Daltons
	if(m_lErrorType & T_FRAGMENT_PPM) fPpmUser = (float)(m_fErr*1e6/200.0); //convert back to true ppm parameter
	else fPpmUser = m_fErr;

	//Check out twice as far for half the score - this is equivalent to checking the next bin
	fPpmUser2x = (float)(fPpmUser*2.0);

	int i=1;
	int iBin;
	
	//perform a single pass through each array.
	//check every point in m_pfSeq, but don't revisit positions in m_vmiType
	for (int a = 0; m_pfSeq[a] != 0; a++) {

    float fIon = (float) m_pfSeq[a];

		while(fIon > m_vmiType[m_lId][i].m_fI){
			i++;
			if(i==m_vmiType[m_lId].size()){
				i--;
				break;
			}
		}
		if(fIon > m_vmiType[m_lId][i].m_fI+1.0) break;

		//Use different calculation based on unit type
		if( (m_vmiType[m_lId][i].m_fI-fIon) < (fIon-m_vmiType[m_lId][i-1].m_fI) ){
			if(m_lErrorType & T_FRAGMENT_PPM)	fPpm = (float)(-(fIon-m_vmiType[m_lId][i].m_fI)/m_vmiType[m_lId][i].m_fI*1e6);
			else fPpm = m_vmiType[m_lId][i].m_fI - fIon;
			iBin=i;
		} else {
			if(m_lErrorType & T_FRAGMENT_PPM) fPpm = (float)((fIon-m_vmiType[m_lId][i-1].m_fI)/m_vmiType[m_lId][i-1].m_fI*1e6);
			else fPpm = fIon - m_vmiType[m_lId][i-1].m_fI;
			iBin=i-1;
		}

		//Check within tolerance
		if(fPpm<fPpmUser){
			iTmp = m_vmiType[m_lId][iBin].m_lM;
      if ((int)m_miUsed.get(iBin) < iTmp) {
        dScore += (double)(iTmp - (int)m_miUsed.get(iBin));
        m_miUsed[iBin] = (float)iTmp;
        lCount++;
      }
			continue;
    }

		//Check at twice the tolerance
		if(fPpm<fPpmUser2x){
			iTmp = m_vmiType[m_lId][iBin].m_lM / 2;
      if ((int)m_miUsed.get(iBin) < iTmp) {
        dScore += (double)(iTmp - (int)m_miUsed.get(iBin));
        m_miUsed[iBin] = (float)iTmp;
      }
    }

	}

	dScore/=1000.0; //this returns dScore to something like its original incarnation.
	*_v = lCount;    
  return (dScore);
}
/*
* MH: The add_N functions are overloaded to support calculation of high res
* mass fragments in the context of the scoring algorithm. m_pfSeq holds
* the accurate mass fragments.
*/
bool mscore_hrk::add_A(const unsigned long _t,const long _c)
{
	unsigned long a = 0;

  //get the conversion factor between a straight sequence mass and an a-ion
	double dValue = m_pSeqUtilFrag->m_dA;
 
	//deal with protein N-terminus
	if(m_bIsN) dValue += m_pSeqUtilFrag->m_fNT;		

	//deal with non-hydrolytic cleavage
	dValue += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);
	if(m_Term.m_lN)	dValue += m_pSeqUtilFrag->m_pdAaMod['['];
	dValue += m_pSeqUtilFrag->m_pdAaFullMod['['];
	unsigned long lValue = 0;

	//calculate the conversion factor between an m/z value and its integer value
	//as referenced in m_vsmapMI
	char cValue = '\0';
	float *pfScore = m_pSeqUtilFrag->m_pfAScore;
	unsigned long lCount = 0;

	//from N- to C-terminus, calcuate fragment ion m/z values and store the results
	//look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a < m_lSeqLength)	{
		cValue = m_pSeq[a];
		dValue += m_pSeqUtilFrag->getAaMass(cValue, tPos+a);
		lValue = mconvert(dValue, _c);
		m_plSeq[lCount] = lValue;

		//MH: tandem uses m_pfSeq differently than K-score.
		//m_pfSeq[lCount] = dValue/(double)_c+m_pSeqUtilFrag->m_dProton;
		m_pfSeq[lCount] = (float)((dValue+_c*m_pSeqUtilFrag->m_dProton)/(double)_c);
		lCount++;
		a++;
	}
 
	//set the next integer mass value to 0: this marks the end of the array 
	m_plSeq[lCount] = 0;
	m_pfSeq[lCount] = 0.0;
	return true;
}

bool mscore_hrk::add_B(const unsigned long _t,const long _c)
{
	unsigned long a = 0;
	
	//get the conversion factor between a straight sequence mass and a b-ion
	double dValue = m_pSeqUtilFrag->m_dB;

	//deal with protein N-terminus
	if(m_bIsN) dValue += m_pSeqUtilFrag->m_fNT;		

	//deal with non-hydrolytic cleavage
	dValue += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);
	if(m_Term.m_lN)	dValue += m_pSeqUtilFrag->m_pdAaMod['['];
	dValue += m_pSeqUtilFrag->m_pdAaFullMod['['];
	unsigned long lValue = 0;

	//calculate the conversion factor between an m/z value and its integer value
	//as referenced in m_vsmapMI
	char cValue = '\0';
	long lCount = 0;
	float *pfScore = m_pSeqUtilFrag->m_pfBScore;
	float *pfScorePlus = m_pSeqUtilFrag->m_pfYScore;

	//from N- to C-terminus, calcuate fragment ion m/z values and store the results
	//look up appropriate scores from m_pSeqUtilFrag->m_pfBScore
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a < m_lSeqLength-1)	{
		cValue = m_pSeq[a];
		dValue += m_pSeqUtilFrag->getAaMass(cValue, tPos+a);
		lValue = mconvert(dValue, _c);
		m_plSeq[lCount] = lValue;

		//MH: tandem uses m_pfSeq differently than K-score.
		m_pfSeq[lCount] = (float)((dValue+_c*m_pSeqUtilFrag->m_dProton)/(double)_c);
		lCount++;
		a++;
	}
	m_plSeq[lCount] = 0;
	m_pfSeq[lCount] = 0.0;
	return true;
}

bool mscore_hrk::add_C(const unsigned long _t,const long _c)
{
	unsigned long a = 0;

	//get the conversion factor between a straight sequence mass and a b-ion
	double dValue = m_pSeqUtilFrag->m_dC;

	//deal with protein N-terminus
	if(m_bIsN) dValue += m_pSeqUtilFrag->m_fNT;		

	//deal with non-hydrolytic cleavage
	dValue += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);
	if(m_Term.m_lN)	dValue += m_pSeqUtilFrag->m_pdAaMod['['];
	dValue += m_pSeqUtilFrag->m_pdAaFullMod['['];
	unsigned long lValue = 0;

	//calculate the conversion factor between an m/z value and its integer value
	//as referenced in m_vsmapMI
	char cValue = '\0';
	long lCount = 0;
	float *pfScore = m_pSeqUtilFrag->m_pfBScore;
	float *pfScorePlus = m_pSeqUtilFrag->m_pfYScore;

	//from N- to C-terminus, calcuate fragment ion m/z values and store the results
	//look up appropriate scores from m_pSeqUtilFrag->m_pfBScore
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a < m_lSeqLength-2)	{
		cValue = m_pSeq[a];
		dValue += m_pSeqUtilFrag->getAaMass(cValue, tPos+a);
		lValue = mconvert(dValue, _c);
		m_plSeq[lCount] = lValue;

		//MH: tandem uses m_pfSeq differently than K-score.
		m_pfSeq[lCount] = (float)((dValue+_c*m_pSeqUtilFrag->m_dProton)/(double)_c);
		lCount++;
		a++;
	}
	m_plSeq[lCount] = 0;
	m_pfSeq[lCount] = 0.0;
	return true;
}

bool mscore_hrk::add_X(const unsigned long _t,const long _c)
{
	long a = m_lSeqLength - 1;

	//get the conversion factor between a straight sequence mass and an x-ion
	double dValue = m_pSeqUtilFrag->m_dX;

	//deal with non-hydrolytic cleavage
	dValue += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	dValue += m_pSeqUtilFrag->m_pdAaMod[']'];
	dValue += m_pSeqUtilFrag->m_pdAaFullMod[']'];

	//deal with protein C-teminus
	if(m_bIsC) dValue += m_pSeqUtilFrag->m_fCT;		
	unsigned long lValue = 0;

	//calculate the conversion factor between an m/z value and its integer value
	//as referenced in m_vsmapMI
	char cValue = '\0';
	unsigned long lCount = 0;
	float fSub = 0.0;
	float *pfScore = m_pSeqUtilFrag->m_pfXScore;

	//from C- to N-terminus, calcuate fragment ion m/z values and store the results
	//look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a > 0)	{
		cValue = m_pSeq[a];
		dValue += m_pSeqUtilFrag->getAaMass(cValue, tPos+a);
		lValue = mconvert(dValue, _c);
		m_plSeq[lCount] = lValue;

		//MH: tandem uses m_pfSeq differently than K-score.
		m_pfSeq[lCount] = (float)((dValue+_c*m_pSeqUtilFrag->m_dProton)/(double)_c);
		lCount++;
		a--;
	}

	//set the next integer mass value to 0: this marks the end of the array 
	m_plSeq[lCount] = 0;
	m_pfSeq[lCount] = 0.0;
	return true;
}

bool mscore_hrk::add_Y(const unsigned long _t,const long _c)
{
	long a = m_lSeqLength - 1;

	//get the conversion factor between a straight sequence mass and a y-ion
	double dValue = m_pSeqUtilFrag->m_dY;
	unsigned long lValue = 0;

	//deal with non-hydrolytic cleavage
	dValue += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	dValue += m_pSeqUtilFrag->m_pdAaMod[']'];
	dValue += m_pSeqUtilFrag->m_pdAaFullMod[']'];

	//deal with protein C-teminus
	if(m_bIsC) dValue +=  m_pSeqUtilFrag->m_fCT;		
	char cValue = '\0';
	unsigned long lCount = 0;
	float fSub = 0.0;
	float *pfScore = m_pSeqUtilFrag->m_pfYScore;
	float *pfScoreMinus = m_pSeqUtilFrag->m_pfBScore;

	//from C- to N-terminus, calcuate fragment ion m/z values and store the results
	//look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a > 0)	{
		cValue = m_pSeq[a];
		dValue += m_pSeqUtilFrag->getAaMass(cValue, tPos+a);
		lValue = mconvert(dValue, _c);
		if(_t == 0)	{
			if(a < 5)	{
				m_plSeq[lCount] = lValue;

				//MH: tandem uses m_pfSeq differently than K-score.
				m_pfSeq[lCount] = (float)((dValue+_c*m_pSeqUtilFrag->m_dProton)/(double)_c);
				lCount++;
			}
		}	else	{
			m_plSeq[lCount] = lValue;
			//MH: tandem uses m_pfSeq differently than K-score.
			m_pfSeq[lCount] = (float)((dValue+_c*m_pSeqUtilFrag->m_dProton)/(double)_c);
			lCount++;
		}
		a--;
	}

	//set the next integer mass value to 0: this marks the end of the array 
	m_plSeq[lCount] = 0;
	m_pfSeq[lCount] = 0.0;
	return true;
}


bool mscore_hrk::add_Z(const unsigned long _t,const long _c)
{
	long a = m_lSeqLength - 1;

	//get the conversion factor between a straight sequence mass and a y-ion
	double dValue = m_pSeqUtilFrag->m_dZ;
	unsigned long lValue = 0;

	//deal with non-hydrolytic cleavage
	dValue += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	dValue += m_pSeqUtilFrag->m_pdAaMod[']'];
	dValue += m_pSeqUtilFrag->m_pdAaFullMod[']'];

	//deal with protein C-teminus
	if(m_bIsC) dValue +=  m_pSeqUtilFrag->m_fCT;		
	char cValue = '\0';
	unsigned long lCount = 0;
	float fSub = 0.0;
	float *pfScore = m_pSeqUtilFrag->m_pfYScore;
	float *pfScoreMinus = m_pSeqUtilFrag->m_pfBScore;

	//from C- to N-terminus, calcuate fragment ion m/z values and store the results
	//look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a > 0)	{
		cValue = m_pSeq[a];
		dValue += m_pSeqUtilFrag->getAaMass(cValue, tPos+a);
		lValue = mconvert(dValue, _c);
		m_plSeq[lCount] = lValue;

		//MH: tandem uses m_pfSeq differently than K-score.
		m_pfSeq[lCount] = (float)((dValue+_c*m_pSeqUtilFrag->m_dProton)/(double)_c);
		lCount++;
		a--;
	}

	//set the next integer mass value to 0: this marks the end of the array 
	m_plSeq[lCount] = 0;
	m_pfSeq[lCount] = 0;
	return true;
}
