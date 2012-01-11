/*
 Copyright (C) 2003-2004 Ronald C Beavis, all rights reserved
 X! tandem 
 This software is a component of the X! proteomics software
 development project

Use of this software governed by the Artistic license, as reproduced here:

The Artistic License for all X! software, binaries and documentation

Preamble
The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some 
semblance of artistic control over the development of the package, 
while giving the users of the package the right to use and distribute 
the Package in a more-or-less customary fashion, plus the right to 
make reasonable modifications. 

Definitions
"Package" refers to the collection of files distributed by the Copyright 
	Holder, and derivatives of that collection of files created through 
	textual modification. 

"Standard Version" refers to such a Package if it has not been modified, 
	or has been modified in accordance with the wishes of the Copyright 
	Holder as specified below. 

"Copyright Holder" is whoever is named in the copyright or copyrights 
	for the package. 

"You" is you, if you're thinking about copying or distributing this Package. 

"Reasonable copying fee" is whatever you can justify on the basis of 
	media cost, duplication charges, time of people involved, and so on. 
	(You will not be required to justify it to the Copyright Holder, but 
	only to the computing community at large as a market that must bear 
	the fee.) 

"Freely Available" means that no fee is charged for the item itself, 
	though there may be fees involved in handling the item. It also means 
	that recipients of the item may redistribute it under the same
	conditions they received it. 

1. You may make and give away verbatim copies of the source form of the 
Standard Version of this Package without restriction, provided that 
you duplicate all of the original copyright notices and associated 
disclaimers. 

2. You may apply bug fixes, portability fixes and other modifications 
derived from the Public Domain or from the Copyright Holder. A 
Package modified in such a way shall still be considered the Standard 
Version. 

3. You may otherwise modify your copy of this Package in any way, provided 
that you insert a prominent notice in each changed file stating how and 
when you changed that file, and provided that you do at least ONE of the 
following: 

a.	place your modifications in the Public Domain or otherwise make them 
	Freely Available, such as by posting said modifications to Usenet 
	or an equivalent medium, or placing the modifications on a major 
	archive site such as uunet.uu.net, or by allowing the Copyright Holder 
	to include your modifications in the Standard Version of the Package. 
b.	use the modified Package only within your corporation or organization. 
c.	rename any non-standard executables so the names do not conflict 
	with standard executables, which must also be provided, and provide 
	a separate manual page for each non-standard executable that clearly 
	documents how it differs from the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

4. You may distribute the programs of this Package in object code or 
executable form, provided that you do at least ONE of the following: 

a.	distribute a Standard Version of the executables and library files, 
	together with instructions (in the manual page or equivalent) on 
	where to get the Standard Version. 
b.	accompany the distribution with the machine-readable source of the 
	Package with your modifications. 
c.	give non-standard executables non-standard names, and clearly 
	document the differences in manual pages (or equivalent), together 
	with instructions on where to get the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

5. You may charge a reasonable copying fee for any distribution of 
this Package. You may charge any fee you choose for support of 
this Package. You may not charge a fee for this Package itself. 
However, you may distribute this Package in aggregate with other 
(possibly commercial) programs as part of a larger (possibly 
commercial) software distribution provided that you do not a
dvertise this Package as a product of your own. You may embed this 
Package's interpreter within an executable of yours (by linking); 
this shall be construed as a mere form of aggregation, provided that 
the complete Standard Version of the interpreter is so embedded. 

6. The scripts and library files supplied as input to or produced as 
output from the programs of this Package do not automatically fall 
under the copyright of this Package, but belong to whomever generated 
them, and may be sold commercially, and may be aggregated with this 
Package. If such scripts or library files are aggregated with this 
Package via the so-called "undump" or "unexec" methods of producing 
a binary executable image, then distribution of such an image shall 
neither be construed as a distribution of this Package nor shall it 
fall under the restrictions of Paragraphs 3 and 4, provided that you 
do not represent such an executable image as a Standard Version of 
this Package. 

7. C subroutines (or comparably compiled subroutines in other languages) 
supplied by you and linked into this Package in order to emulate 
subroutines and variables of the language defined by this Package 
shall not be considered part of this Package, but are the equivalent 
of input as in Paragraph 6, provided these subroutines do not change 
the language in any way that would cause it to fail the regression 
tests for the language. 

8. Aggregation of this Package with a commercial distribution is always 
permitted provided that the use of this Package is embedded; that is, 
when no overt attempt is made to make this Package's interfaces visible 
to the end user of the commercial distribution. Such use shall not be 
construed as a distribution of this Package. 

9. The name of the Copyright Holder may not be used to endorse or promote 
products derived from this software without specific prior written permission. 

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE. 

The End 
*/

// File version: 2005-01-23

#include "stdafx.h"
#include "msequence.h"
#include "mspectrum.h"
#include "msequtilities.h"
#include "mscore_tandem.h"

// Factory instance, registers itself with the mscoremanager.
static mscorefactory_tandem factory;
	
mscorefactory_tandem::mscorefactory_tandem()
{
	mscoremanager::register_factory("tandem", this);
}

mplugin* mscorefactory_tandem::create_plugin()
{
	return new mscore_tandem();
}

mscore_tandem::mscore_tandem(void)
{
	m_dScale = 4.0;
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

mscore_tandem::~mscore_tandem(void)
{
	if(m_pFactorial != NULL)
		delete m_pFactorial;
}

/*
 * add_mi does the work necessary to set up an mspectrum object for modeling. 
 *   - an entry in the m_State object is made for the parent ion M+H
 * once an mspectrum has been added, the original mspectrum is no longer
 * needed for modeling, as all of the work associated with a spectrum
 * is only done once, prior to modeling sequences.
 */
bool mscore_tandem::add_mi(mspectrum &_s)
{
	if (!mscore::add_mi(_s))
		return false;

	if (m_vmiType.size() == 0)
		m_vmiType.reserve((long)m_vSpec.size());
/*
 * use blur to improve accuracy at the edges of the fragment ion m/z error range
 */
	blur(_s.m_vMI);

	return true;
}

/*
 * blur takes an SMap object and adds values to it adjacent to the 
 * values in the existing map. this process makes the conversion from
 * floating point m/z values to integer SMap index values more accurate.
 * the half-width of the distribution around each initial value is set
 * by the m_fWidth value. For example, 
 *
 * if there is an intensity value I at index value N, 
 * then with a width of 2, the following new values are created,
 * I(N-2),I(N-1),I(N+1),I(N+2)
 *
 * the value for the width has an impact on the performance of the
 * protein modeling session: the wider the width, the more memory that is required to hold
 * the SMaps for the spectra. as that memory expands, the ability of
 * the various processor caches to hold the SMaps will change. keeping
 * the SMaps in the fastest cache speeds up the calculation considerably.
 */
inline bool mscore_tandem::blur(vector<mi> &_s)
{
	vmiType vType;
	MIType uType;
	uType.m_fI = 0.0;
	uType.m_lM = 0;
	long lValue = 0;
	long a = 0;
	long w = (long)(0.1+m_fWidth);
/*
 * fConvert and fFactor are only used if the m_lErrorType uses ppm for fragment ion mass errors
 * if ppm is used, the width at m/z = 200.0 is taken as the base width for blurring & the
 * width is scaled by the measured m/z. 
 * NOTE: the m_fErr value used in the ppm case is: (the error in ppm) x 200
 */
	float fConvert = m_fErr/m_fWidth;
	const float fFactor = (float)200.0/fConvert;
	const size_t tSize = _s.size();
	size_t tCount = 0;
	vType.reserve(tSize*3);
/*
 * add additional values based on the m_fWidth setting
 */
	while(tCount < tSize)	{
		if(_s[tCount].m_fI > 0.5)	{
			lValue = (long)(_s[tCount].m_fM/fConvert);
			a = -1*w;
			if(m_lErrorType & T_FRAGMENT_PPM)	{
				a = (long)((float)a * (float)lValue/fFactor - 0.5);
				if(a > -1*w)
					a = -1*w;
			}
			while(a <= w)	{
				if(uType.m_lM == lValue + a)	{
					if(uType.m_fI < _s[tCount].m_fI)	{
						vType.back().m_fI = _s[tCount].m_fI;
					}
				}
				else	{
						uType.m_lM = lValue + a;
						uType.m_fI = _s[tCount].m_fI;
						vType.push_back(uType);
				}
				a++;
			}
		}
		tCount++;
	}
	m_vmiType.push_back(vType);
	return true;
}

/*
 * hfactor returns a factor applied to the score to produce the
 * hyper score, given a number of ions matched.
 */
double mscore_tandem::hfactor(long _l) {
	return m_pFactorial[_l];
}
/*
 * hconvert is use to convert a hyper score to the value used in a score
 * histogram, since the actual scoring distribution may not lend itself
 * to statistical analysis through use of a histogram.
 */
float mscore_tandem::hconvert(float _f) {
	if(_f <= 0.0)
		return 0.0;
	return (float)(m_dScale*log10(_f));
}

bool mscore_tandem::clear()
{
	m_vmiType.clear();
	return true;
}

/*
 * dot is the fundamental logic for scoring a peptide with a mass spectrum.
 * the mass spectrum is determined by the value of m_lId, which is its index
 * number in the m_vsmapMI vector. the sequence is represented by the values
 * that are currently held in m_plSeq (integer masses) and m_pfSeq (scoring
 * weight factors). 
 * convolution scores are the sum of the products of the spectrum intensities
 * and the scoring weight factors for a spectrum.
 * hyper scores are the product of the products of the spectrum intensities
 * and the scoring weight factors for a spectrum.
 * the number of predicted ions that correspond to non-zero intensities in
 * the spectrum are tallied and returned.
 */

double mscore_tandem::dot(unsigned long *_v)
{
	float fScore = 0.0;
	float fValue0 = 0.0;
	unsigned long a = 0;
	unsigned long lCount = 0;
	long lType = 0;
	size_t b = 0;
	vector<MIType>::iterator itType = m_vmiType[m_lId].begin();
	vector<MIType>::const_iterator itEnd = m_vmiType[m_lId].end();
	// tType and tTypeSize were added in 2006.09.01 to correct a problem
	// created by VC++ 2005. This new version uses a strict bounds checking
	// style for STL iterators, that cause a run time error if an iterator
	// longer than .end() is produced by incrementing the iterator.
	size_t tType = 0;
	size_t tTypeSize = m_vmiType[m_lId].size();
	const size_t tStep = 5;
	while(m_plSeq[a] != 0 && itType != itEnd)	{
		lType = 0;
		if(itType->m_lM < m_plSeq[a])	{
			lType = 1;
			// Generally lots more spectrum peaks than sequence peaks.  Trying
			// large steps first helps reduce performance hit for this operation.
			// This is were the iterator bounds checking failed in VC++ 2005.
			// By checking the size first, the iterator is not evaluated and
			// does not produce the failure.
			while (tType+tStep < tTypeSize && (itType+tStep)->m_lM < m_plSeq[a]) {
				itType += tStep;
				tType += tStep;
			}
			do {
				itType++;
				tType++;
			} while(itType != itEnd && itType->m_lM < m_plSeq[a]);
		}
		else if(itType->m_lM > m_plSeq[a])	{
			do {
				a++;
			} while(itType->m_lM > m_plSeq[a] && m_plSeq[a] != 0);
		}
		if(m_plSeq[a] == 0 || itType == itEnd)	{
			break;
		}
		if(itType->m_lM == m_plSeq[a])	{
			fValue0 = itType->m_fI * m_pfSeq[a];
			if(fValue0 > 0.0)	{
				lCount++;
				fScore += fValue0;
			}
		}
		if(lType)	{
			a++;
		}
		else	{
			itType++;
			tType++;
		}
	}
	*_v = lCount;	
	return (fScore);
}
