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

#include "mscore_kgpu.h"
#include "cuda.h"
#include <iostream>
#include <sstream>

// Factory instance, registers itself with the mscoremanager.
static mscorefactory_kgpu factory;
static bool first_access = true;

// format a byte count in a human readable manner
std::string prettymem(size_t bytes) {
  std::ostringstream result;
  result << bytes << " bytes";
  if (bytes >= (1024*1024)) {
	result << " (" << bytes/(1024*1024) <<" MB)";
  }
  else if (bytes >= 1024 ) {
	result << " (" <<bytes/1024 <<" KB)";
  }
  return result.str();
}

// logging stuff, stolen from FastPass project
enum LogType {logall, debug, info, warn, error, lognone};
enum EventType {noevent, init, initok, hostalloc, hostrelease, hostreleaseok, hostallocok, devicealloc, devicerelease, devicereleaseok, deviceallocok, parse, parseok, copye, copyok, kernel, kernelok, file, fileok};

// ug, better design is wanted here
static LogType logLevel=debug;

#define logevent(LOGTYPE, EVENT, LOGOUTPUT)	\
  if (LOGTYPE >= logLevel) { \
  std::cerr << #LOGTYPE << "(" << #EVENT << "): " << LOGOUTPUT ; \
  }

// end logging stuff


mscorefactory_kgpu::mscorefactory_kgpu()
{
    m_initialFreeMemory = 0;
    mscoremanager::register_factory("kgpu-score", this);
}

mplugin* mscorefactory_kgpu::create_plugin()
{
    if (first_access) { // first call, from master thread
        first_access = false; // just once, thanks
        // initialize GPU
        int deviceCount = 0;
        int bestDeviceIndex = 0;
        double highestDeviceScore = 0.0;  // device score = clockRate * multiProcessorCount
        cudaDeviceProp deviceProp;

        cuInit(0);

        // loop through all devices and choose which one to use based on
        // highest device score

        cudaError_t err = cudaGetDeviceCount(&deviceCount);
        if (err != cudaSuccess) {
             logevent(error, init, cudaGetErrorString(err) << endl);
        }
        logevent(info, init, "GPU devices found: " << deviceCount << endl);

        for (int curDeviceNumber=0; curDeviceNumber < deviceCount; curDeviceNumber++) {
            logevent(info, init, "testing device " << curDeviceNumber << endl);
            memset(&deviceProp, 0, sizeof(deviceProp));
            if (cudaSuccess == cudaGetDeviceProperties(&deviceProp, curDeviceNumber)) {
                double curDeviceScore = deviceProp.multiProcessorCount * deviceProp.clockRate;
                if (curDeviceScore > highestDeviceScore) {
                    highestDeviceScore = curDeviceScore;
                    bestDeviceIndex = curDeviceNumber;
                }
            }
            else {
                logevent(error, init, "unable to access GPU properites, using non-GPU k-score" << endl);
                return new mscore_k();
            }
        }

        memset(&deviceProp, 0, sizeof(deviceProp));
        if (cudaSuccess == cudaGetDeviceProperties(&deviceProp, bestDeviceIndex)) {
            logevent(info, init, "Using GPU device " << bestDeviceIndex << ": \"" << deviceProp.name << "\", "<< prettymem(deviceProp.totalGlobalMem) << ", " << deviceProp.multiProcessorCount << " multiprocessors, " << deviceProp.clockRate/1000 << " MHz" << endl);
            size_t free_mem, total_mem;
            cudaMemGetInfo(&free_mem, &total_mem);
            m_initialFreeMemory = free_mem;
            logevent(info, init, "device initial free memory: " << prettymem(m_initialFreeMemory) << endl);
        }
        else {
            logevent(error, init, "unable to access device properites, using non-GPU k-score" << endl);
            return new mscore_k();
        }      

        logevent(info, init, "resetting GPU" << endl);
        if (cudaSuccess != cudaDeviceReset() ) {
            logevent(error, init, "error resetting GPU" << endl);
        }
        logevent(info, initok, "done resetting GPU" << endl);

        if (cudaSuccess != cudaSetDevice(bestDeviceIndex)) {
            logevent(error, init, "Error - cannot set GPU, using non-GPU k-score" << bestDeviceIndex << endl);
            return new mscore_k();
        }
        else {
            logevent(info, initok, "GPU " << bestDeviceIndex << " initialized." << endl);
        }
        size_t free, total;
        free = total = 0;
        cudaMemGetInfo(&free, &total);
        m_initialFreeMemory = total;
        logevent(info, init, "Device memory: " << prettymem(free) << " free / " << prettymem(m_initialFreeMemory) << " total" << endl);

    }
    if (!m_initialFreeMemory) {
        return new mscore_k();
    }
    return new mscore_kgpu();
}

mscore_kgpu::mscore_kgpu(void)
{
    m_dScale = 0.05;
    m_maxEnd = 0;
    m_dIsotopeCorrection = 1.0;
}

mscore_kgpu::~mscore_kgpu(void)
{
}

bool mscore_kgpu::clear()
{
    m_vmiType.clear();
    return true;
}

/*
 * allows score object to issue warnings, or set variable based on xml.
 */
bool mscore_kgpu::load_param(XmlParameter &_x)
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
bool mscore_kgpu::precondition(mspectrum &_s)
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
void mscore_kgpu::prescore(const size_t _i)
{
    mscore::prescore(_i);

    // Initialize of clear the used intensity look-up.
    if (m_miUsed.m_pfI == NULL)
        m_miUsed.init(0, m_maxEnd);
    else
        m_miUsed.clear();
}

/*
 * add_mi does the work necessary to set up an mspectrum object for modeling. 
 *   - an entry in the m_State object is made for the parent ion M+H
 * once an mspectrum has been added, the original mspectrum is no longer
 * needed for modeling, as all of the work associated with a spectrum
 * is only done once, prior to modeling sequences.
 */
bool mscore_kgpu::add_mi(mspectrum &_s)
{
    if (!mscore::add_mi(_s))
        return false;

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
 * mconvert converts from mass and charge to integer ion value
 * for mi vector.
 */
unsigned long mscore_kgpu::mconvert(double _m, const long _c)
{
    const double fZ = (double)_c;

    double dMass = (fZ*m_pSeqUtilFrag->m_dProton + _m)/fZ;
    return imass(dMass);
}

/*
 * sfactor returns a factor applied to the final convolution score.
 */
double mscore_kgpu::sfactor()
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
void mscore_kgpu::report_score(char* buffer, float hyperscore)
{
    sprintf(buffer, "%d",(int) (hyperscore + 0.5));
}

/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses).
 */
double mscore_kgpu::dot(unsigned long *_v)
{
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
						if (m_miUsed[iIon] < fTmp) {
								dScore += fTmp - m_miUsed[iIon];
								m_miUsed[iIon] = fTmp;
								lCount++;
						}

						itGreater = itType + 1;
				}

				if (itType != itBegin && iIon - 1 == itType[-1].m_lM) {
						fTmp = ((float)0.5)*itType[-1].m_fI;
						if (m_miUsed[iIon-1] < fTmp) {
								dScore += fTmp -m_miUsed[iIon-1];
								m_miUsed[iIon-1] = fTmp;
						}
				}

				if (itGreater != itEnd && iIon + 1 == itGreater->m_lM) {
						fTmp = ((float)0.5)*itGreater->m_fI;
						if (m_miUsed[iIon+1] < fTmp) {
								dScore += fTmp - m_miUsed[iIon+1];
								m_miUsed[iIon+1] = fTmp;
						}
				}
		}

    *_v = lCount;    
    return (dScore);
}

