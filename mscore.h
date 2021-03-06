// -*- mode: c++ -*-

/*
 Copyright (C) 2003 Ronald C Beavis, all rights reserved
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

#ifndef MSCORE_H
#define MSCORE_H
/*
	The pluggable scoring system was devised and implemented by Brendan Maclean. It allows
	developers to simply add new scoring systems to X! Tandem for any purpose. Unless you have
	some specific need, it is best to leave PLUGGABLE_SCORING defined.
*/
//#define PLUGGABLE_SCORING
// File version: 2004-02-01
// File version: 2004-03-01

/* 
   changes for MPI and Hadoop bpratt Insilicos LLC after:
   Modified 2007 Robert D Bjornson for X!!Tandem, Parallel MPI version.
*/

#include <set>

#include "mscorestate.h"
#include "mscorepam.h"
#include "mplugin.h"
#include "msequtilities.h"

class XmlParameter;

class mprocess; // make sure the load_next() declaration understands that mprocess is a class to avoid problems with casts, otherwise.  


class mspectrumindex
{
public:
	mspectrumindex() { }
	virtual ~mspectrumindex() { }
	
	float m_fM; // the M+H + error for an mspectrum
	unsigned long m_tA; // the index number for an mspectrum, in the m_vSpectra vector
/*
 * override the less than operator, so that an mspectrumdetails can be easily compared to a
 * float M+H value
 */
	bool operator<(const mspectrumindex &rhs) const
	{	
		return m_fM < rhs.m_fM; 
	}
	bool operator<=(const mspectrumindex &rhs) const
	{	
		return m_fM <= rhs.m_fM; 
	}
	bool operator>=(const mspectrumindex &rhs) const
	{	
		return m_fM >= rhs.m_fM; 
	}
/*
 * override the greater than operator, so that an mspectrumdetails can be easily compared to a
 * float M+H value
 */
	bool operator>(const mspectrumindex &rhs) const
	{	
		return m_fM > rhs.m_fM; 
	}
/*
 * override the equivalence operator, so that an mspectrumdetails can be easily compared to a
 * float M+H value
 */
	bool operator ==(const mspectrumindex &rhs) const
	{	
		return (rhs.m_fM == m_fM); 
	}
/*
 * simple copy operation, using the = operator
 */
	mspectrumindex& operator=(const mspectrumindex &rhs)
	{
		m_fM = rhs.m_fM;
		m_tA = rhs.m_tA;
		return *this;
	}
#ifdef HAVE_MULTINODE_TANDEM  // support for Hadoop and/or MPI?
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
		ar & m_fM;
		ar & m_tA;
    }
#endif
};

/*
 * mspectrumdetails is a specialty class used by mscore to check whether a particular mspectrum
 * has a parent ion M+H within error of a given peptide sequence (with modifications if appropriate)
 */
class mspectrumdetails
{
public:
	mspectrumdetails() { }
	virtual ~mspectrumdetails() { }
	
	float m_fU; // the M+H + error for an mspectrum
	float m_fL; // the M+H - error for an mspectrum
	long m_lA; // the index number for an mspectrum, in the m_vSpectra vector
#ifdef HAVE_MULTINODE_TANDEM  // support for Hadoop and/or MPI?
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
		ar & m_fU;
		ar & m_fL;
		ar & m_lA;
    }
#endif
  /*
 * override the less than operator, so that an mspectrumdetails can be easily compared to a
 * float M+H value
 */
	bool operator<(const float &rhs)
	{	
		return m_fL < rhs; 
	}
/*
 * override the greater than operator, so that an mspectrumdetails can be easily compared to a
 * float M+H value
 */
	bool operator>(const float &rhs)
	{	
		return m_fU > rhs; 
	}
/*
 * override the equivalence operator, so that an mspectrumdetails can be easily compared to a
 * float M+H value
 */
	bool operator ==(const float &rhs)
	{	
		return (rhs >= m_fL && rhs <= m_fU); 
	}
/*
 * simple copy operation, using the = operator
 */
	mspectrumdetails& operator=(const mspectrumdetails &rhs)
	{
		m_fU = rhs.m_fU;
		m_fL = rhs.m_fL;
		m_lA = rhs.m_lA;
		return *this;
	}
};

class PermuteState
{
public:
	PermuteState() { m_pSeq = new char[256]; m_pPerm = new char[256]; m_lSize = 255;
					 m_tPos = 0; m_tEnd = 0;*m_pSeq=0;*m_pPerm=0;m_bRev=false;} // bpratt 4-5-11
	~PermuteState() { delete[] m_pSeq; delete[] m_pPerm;}
	size_t m_tPos;
	size_t m_tEnd;
	char *m_pSeq;
	char *m_pPerm;
	unsigned long m_lSize;
	bool m_bRev;
#ifdef HAVE_MULTINODE_TANDEM  // support for Hadoop and/or MPI?
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
      ar & m_tPos;
      ar & m_tEnd;
      ar & m_lSize;
	  size_t n=0;
	  do {
        ar & m_pSeq[n];
	  } while (m_pSeq[n++]);
	  n=0;
	  do {
        ar & m_pPerm[n];
	  } while (m_pPerm[n++]);
      ar & m_bRev;
    }
#endif
};

class MIType
{
#ifdef HAVE_MULTINODE_TANDEM  // support for Hadoop and/or MPI?
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
		ar & m_lM;
		ar & m_fI;
		ar & m_lA;
    }
#endif
public:
	MIType() { }
	virtual ~MIType() { }
	
	unsigned long m_lM; // the M+H + error for an mspectrum
	float m_fI; // the M+H - error for an mspectrum
	long m_lA; // the index number for an mspectrum, in the m_vSpectra vector
/*
 * simple copy operation, using the = operator
 */
	MIType& operator=(const MIType &rhs)
	{
		m_fI = rhs.m_fI;
		m_lM = rhs.m_lM;
		return *this;
	}
};

typedef vector<MIType> vmiType;
class vectorvmiType : public vector<vmiType> {
#ifdef HAVE_MULTINODE_TANDEM  // autoserialize fails for vector of vectors
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
		size_t s = this->size();
		size_t initial_s = s;
		ar & s;
		if (s > initial_s) { // loading
			this->resize(s);
		}
		while (s--) {
			ar & (*this)[s];
		}
    }
#endif
};

class mspec
{
#ifdef HAVE_MULTINODE_TANDEM
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
		ar & m_fMH;
		ar & m_fZ;
    }
#endif
public:
	mspec() { m_fMH = 0.0; m_fZ = 1.0;}
	virtual ~mspec() { }
	
	double m_fMH; // the M+H value of an mspectrum
	float m_fZ; // the charge value of an mspectrum
/*
 * simple copy operation, using the = operator
 */
	mspec& operator=(const mspec &rhs)
	{
		m_fZ = rhs.m_fZ;
		m_fMH = rhs.m_fMH;
		return *this;
	}
	mspec& operator=(const mspectrum &rhs)
	{
		m_fZ = rhs.m_fZ;
		m_fMH = (float)rhs.m_dMH;
		return *this;
	}
};
/*
 * mscore is the class that contains most of the logic for comparing one sequence with
 * many tandem mass spectra. mprocess contains a comprehensive example of how to use an mscore
 * class. mscore has been optimized for speed, without any specific modifications to take
 * advantage of processor architectures.
 */
class mscore : public mplugin
{
#ifdef HAVE_MULTINODE_TANDEM
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
    {
#define SERIALIZE(ar,var) {std::string varname(#var);ar & varname;ar & var;}
	SERIALIZE(ar, m_fErr); // error for the fragment ions
	SERIALIZE(ar, m_fHomoError);
	SERIALIZE(ar, m_fHyper); // current hyper score
	SERIALIZE(ar, m_fParentErrMinus); // error for the parent ion M+H (not m/z)
	SERIALIZE(ar, m_fParentErrPlus); // error for the parent ion M+H (not m/z)
	SERIALIZE(ar, m_fMaxMass);
	SERIALIZE(ar, m_fMinMass);
	SERIALIZE(ar, m_lMaxCharge); // current parent ion charge
	SERIALIZE(ar, m_lMaxPeaks); // if > 0, the m_lMaxPeaks most intense peaks will be used
	SERIALIZE(ar, m_dScale); // scale for use in hconvert
    SERIALIZE(ar, m_seqUtil); // class contains variables and constants for calculating sequence and
						     // fragment masses
	SERIALIZE(ar, m_seqUtilAvg); // class contains variables and constants for calculating fragment masses
								// based on average atomic masses
	// omitted SERIALIZE(ar, m_pSeqUtilFrag); // pointer to the msequtilities object to use for fragment ion masses
	// omitted SERIALIZE(ar, m_State); // class stores information about the potential modification state machine
	// omitted SERIALIZE(ar, m_Pam); // class stores information about point mutations state machine
 	// omitted SERIALIZE(ar, m_Sap); // class stores information about single amino acid polymorphisms state machine

	// omitted SERIALIZE(ar, m_Term); // class stores information about potential modification of the N- & C-terminii
	SERIALIZE(ar, m_plCount);// ion count information, indexed using the mscore_type_a enum
	SERIALIZE(ar, m_pfScore);// convolute score information, indexed using the mscore_type_a enum

	SERIALIZE(ar, m_lType); // current ion type - value from mscore_type
	SERIALIZE(ar, m_bUsePam); // true if the peptide will be checked for all possible point mutations
 	SERIALIZE(ar, m_bUseSaps); // true if the peptide will be checked for all known single amino acid polymorphisms
	SERIALIZE(ar, m_bIsC); // true if the current peptide contains the C-terminus of the protein
	SERIALIZE(ar, m_bIsN); // true if the current peptide contains the N-terminus of the protein
	SERIALIZE(ar, m_bIsotopeError); // true if the spectrum mass may be associated with the wrong isotopic peak
	SERIALIZE(ar, m_lMILength); // current length of the mi vector
	SERIALIZE(ar, m_lSeqLength); // current sequence length
	unsigned long former_lsize = m_lSize;
	SERIALIZE(ar, m_lSize); // maximum sequence length - this can be adjusted on the fly
	SERIALIZE(ar, m_lSpectra); // current length of the m_vSpec vector
	SERIALIZE(ar, m_lErrorType); // current ion mass accuracy information - value from mscore_error
	SERIALIZE(ar, m_fScore); // current convolution score
	SERIALIZE(ar, m_dSeqMH); // current sequence M+H - changed from m_fSeqMH to improve accuracy of parent ion mass calculations
	SERIALIZE(ar, m_fWidth); // current half-width of the entry for a single fragment ion in the m_vsmapMI map
	                // this value is used by blur
	if (former_lsize != m_lSize) { // looks like we're loading - initial value is incorrect
		delete[] m_pfSeq;
		delete[] m_plSeq;
		delete[] m_pSeq;
		m_pfSeq = new float[m_lSize];
		m_plSeq = new unsigned long[m_lSize];
		m_pSeq = new char[m_lSize];
	}
	for (size_t a=m_lSize;a--;) {
		ar & m_pSeq[a]; // the current sequence
		ar & m_pfSeq[a]; // residue masses corresponding to the current sequence in daltons
		ar & m_plSeq[a]; // residue masses corresponding to the current sequence, converted into integers
	}
	SERIALIZE(ar, m_lId); // id of the current spectrum
	SERIALIZE(ar, m_tSeqPos); // zero-based absolute position of the current peptide in the protein sequence
	SERIALIZE(ar, m_lDetails);
	SERIALIZE(ar, m_iCharge);
	SERIALIZE(ar, m_bMini);
	SERIALIZE(ar, m_vSpec); // vector of all spectra being considered
	                        // for all spectra being considered
	SERIALIZE(ar, m_vDetails); // vector of mspectrumdetails objects, for looking up parent ion M+H
	                                     // values of mass spectra
	//  huh, std::set doesn't work
	int s = (int)m_sIndex.size();
	ar & s;
	if (s != (int)m_sIndex.size()) { // loading
		while (s--) {
			mspectrumindex itmp;
			ar & itmp;
			m_sIndex.insert(itmp);
		}
	} else {
		set<mspectrumindex>::iterator it = m_sIndex.begin();
		while (it != m_sIndex.end()) {
			mspectrumindex itmp = *it++;
			ar & itmp;
		}
	}
	SERIALIZE(ar, m_psPermute);
    }
#endif
public:
	mscore(void);
	virtual ~mscore(void);
public:
	float m_fErr; // error for the fragment ions
	float m_fHomoError;
	float m_fHyper; // current hyper score
	float m_fParentErrMinus; // error for the parent ion M+H (not m/z)
	float m_fParentErrPlus; // error for the parent ion M+H (not m/z)
	float m_fMaxMass;
	float m_fMinMass;
	long m_lMaxCharge; // current parent ion charge
	unsigned long m_lMaxPeaks; // if > 0, the m_lMaxPeaks most intense peaks will be used
	double m_dScale; // scale for use in hconvert
	msequtilities m_seqUtil; // class contains variables and constants for calculating sequence and
						     // fragment masses
	msequtilities m_seqUtilAvg; // class contains variables and constants for calculating fragment masses
								// based on average atomic masses
	msequtilities* m_pSeqUtilFrag; // pointer to the msequtilities object to use for fragment ion masses
	mscorestate m_State; // class stores information about the potential modification state machine
	mscorepam m_Pam; // class stores information about point mutations state machine
	mscoresap m_Sap; // class stores information about single amino acid polymorphisms state machine
	mscoreterm m_Term; // class stores information about potential modification of the N- & C-terminii
	unsigned long m_plCount[16];// ion count information, indexed using the mscore_type_a enum
	float m_pfScore[16];// convolute score information, indexed using the mscore_type_a enum

	unsigned long m_lType; // current ion type - value from mscore_type

public:

/*
      The following section represents the pluggable scoring API devised and
      implemented by Brendan MacLean. It allows developers to simply add new
      scoring systems to X! Tandem for any purpose.  Please do not modify without
      considering carefully how changes might impact external scoring plug-ins.
*/
      virtual bool load_param(XmlParameter &_x); // allows score object to issue warnings,
                                                                  // or set variables based on xml
      virtual bool precondition(mspectrum &_s); // called before spectrum conditioning 
      virtual bool add_details(mspectrum &_s);
      virtual bool add_mi(mspectrum &_s);
      virtual void prescore(const size_t _i); // called before scoring
      virtual float score(const size_t _i);
/*
      If you have no need for scoring plug-ins other than mscore_tandem, making
      this function non-virtual can yield a 10% perf improvement.
*/
#ifdef PLUGGABLE_SCORING
      virtual unsigned long mconvert(double _m, const long _c); // convert mass to integer ion m/z for mi vector
#else
      unsigned long mconvert(double _m, const long _c); // convert mass to integer ion m/z for mi vector
#endif
      virtual double hfactor(long _l); // hyper scoring factor from number of ions matched
      virtual double sfactor(); // factor applied to final convolution score
      virtual float hconvert(float _h); // convert hyper score to histogram score
      virtual void report_score(char* _buff, float _h); // format hyper score for output
      virtual bool clear();
protected:
      virtual double dot(unsigned long *_v) = 0; // this is where the real scoring happens
/*    End puggable scoring API. */

public:
	unsigned long add_seq(const char *_s,const bool _n,const bool _c,const unsigned long _l,const int _f);
	bool sort_details();
	bool get_aa(vector<maa> &_m,const size_t _a,double &_d);
	virtual bool load_next(const mprocess *_parentProcess);  // made virtual and added arg in support of pluggable scoring bpratt
	bool load_state(void);
#ifdef PLUGGABLE_SCORING
// To mimic SEQUEST-like scoring, a scoring plugin has to modify the theoretical fragment ion series to 
// adjust intensities and add appropriate neutral losses, which is accomplished by overriding this method.
    virtual bool load_seq(const unsigned long _t,const long _c);
#else
	bool load_seq(const unsigned long _t,const long _c);
#endif
	double seq_mh(void); // changed from float to accomodate more accurate parent ion mass calculation (2005.02.01)
	unsigned long set_seq(const char *_s,const bool _n,const bool _c,const unsigned long _l,const int _f);
	bool set_isotope_error(const bool _b);
	unsigned long set_type(const unsigned long _t);
	unsigned long set_error(const unsigned long _t);
	bool set_pam(const bool _b);
	bool set_saps(const bool _b,string &_s);
	float set_parent_error(const float _f,const bool _b);
	float set_homo_error(const float _f);
	float set_fragment_error(const float _f);
	void set_fragment_masstype(masscalc::massType _t);
	bool test_parents(size_t &_t);
	bool set_pos(const size_t _t);
	bool set_mini(const bool _b);
	bool reset_permute();
	bool permute();
public:
	enum	{
		T_Y =	0x01,
		T_B =	0x02,
		T_X =	0x04,
		T_A =	0x08,
		T_C =	0x10,
		T_Z =	0x20,
		T_Y_18 = 0x100,
		T_X_18 = 0x200,
		T_A_18 = 0x400,
		T_B_18 = 0x800,
		T_C_18 = 0x1000,
		T_Y_17 = 0x10000,
		T_X_17 = 0x20000,
		T_A_17 = 0x40000,
		T_B_17 = 0x80000,
		T_C_17 = 0x100000,
	} mscore_type; // enum for referencing information about specific ion types.
	enum	{
		S_Y =	1,
		S_B =	2,
		S_X =	3,
		S_A =	4,
		S_C =	5,
		S_Z =	6,
		S_Y_18 = 7,
		S_X_18 = 8,
		S_A_18 = 9,
		S_B_18 = 10,
		S_Y_17 = 11,
		S_X_17 = 12,
		S_A_17 = 13,
		S_B_17 = 14,
		S_C_17 = 15,
	} mscore_type_a; // enum for referencing information about specific ion types.
	enum	{
		T_PARENT_DALTONS = 0x01,
		T_PARENT_PPM = 0x02,
		T_FRAGMENT_DALTONS = 0x04,
		T_FRAGMENT_PPM =	0x08,
	} mscore_error; // enum for referencing information about ion mass measurement accuracy.

protected:
	bool m_bUsePam; // true if the peptide will be checked for all possible point mutations
	bool m_bUseSaps; // true if the peptide will be checked for all known single amino acid polymorphisms
	bool m_bIsC; // true if the current peptide contains the C-terminus of the protein
	bool m_bIsN; // true if the current peptide contains the N-terminus of the protein
	bool m_bIsotopeError; // true if the spectrum mass may be associated with the wrong isotopic peak
	unsigned long m_lMILength; // current length of the mi vector
	unsigned long m_lSeqLength; // current sequence length
	unsigned long m_lSize; // maximum sequence length - this can be adjusted on the fly
	long m_lSpectra; // current length of the m_vSpec vector
	unsigned long m_lErrorType; // current ion mass accuracy information - value from mscore_error
	float m_fScore; // current convolution score
	double m_dSeqMH; // current sequence M+H - changed from m_fSeqMH to improve accuracy of parent ion mass calculations
	float m_fWidth; // current half-width of the entry for a single fragment ion in the m_vsmapMI map
	                // this value is used by blur
	float *m_pfSeq; // residue masses corresponding to the current sequence in daltons
	unsigned long *m_plSeq; // residue masses corresponding to the current sequence, converted into integers
	char *m_pSeq; // the current sequence
	size_t m_lId; // id of the current spectrum
	size_t m_tSeqPos; // zero-based absolute position of the current peptide in the protein sequence
	long m_lDetails;
	int m_iCharge;
  bool m_bMini;
protected:
	vector<mspec> m_vSpec; // vector of all spectra being considered
	                        // for all spectra being considered
	vector<mspectrumdetails> m_vDetails; // vector of mspectrumdetails objects, for looking up parent ion M+H
	                                     // values of mass spectra
	set<mspectrumindex> m_sIndex;
	PermuteState m_psPermute;

protected:
	virtual bool add_A(const unsigned long _t,const long _c);
	virtual bool add_B(const unsigned long _t,const long _c);
	virtual bool add_C(const unsigned long _t,const long _c);
	virtual bool add_Y(const unsigned long _t,const long _c);
	virtual bool add_X(const unsigned long _t,const long _c);
	virtual bool add_Z(const unsigned long _t,const long _c);
	bool check_parents(void);
	bool load_next_pam(void);
	bool load_next_sap(void);
	bool check_pam_mass();
	bool load_next_term(void);
	bool run_state_machine(void);
};

/*
 * mscoremanager contains static short-cuts for dealing with mscore
 * plug-ins.
 */
class mscoremanager
{
public:
	static const char* TYPE;

	static mscore* create_mscore(XmlParameter &_x);
	static void register_factory(const char* _spec, mpluginfactory* _f);
};

bool lessThanDetails(const mspectrumdetails &_l,const mspectrumdetails &_r);
bool lessThanMI(const mi &_l,const mi &_r);

#endif
