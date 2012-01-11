/*
 MapReduce implementation of X!Tandem
 Copyright (C) 2010 Insilicos LLC All Rights Reserved

 Also includes MPI implementation based on X!!Tandem which is
 Copyright (C) 2007 Robert D Bjornson, all rights reserved
 (The serialization code shared by MapReduce and MPI is based on
 the serialization in the original X!!Tandem work)

 All of which is based on the original X!Tandem project which is
 Copyright (C) 2003-2004 Ronald C Beavis, all rights reserved
 
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

#ifdef HAVE_MULTINODE_TANDEM // support for Hadoop and/or MPI?

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <string>
#include <sys/timeb.h>
#include <ctime>

#include "stdafx.h"
#include "serialize.h"
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "mhistogram.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"
#include "mapreducehelper.h"
#include <unistd.h>


namespace boost {
  namespace serialization {




template<class Archive>
  void serialize(Archive &ar,  msequenceCollection &m, const unsigned int version)
  {
	ar & m.m_tLength;
	ar & m.m_tMax; // maximum number of sequences to be stored in the m_vASequences vector
	ar & m.m_vASequences; // a vector of msequence objects
  }

  template<class Archive>
  void serialize(Archive &ar, XmlParameter &p, const unsigned int version)
    {
	ar & p.m_mapParam; // a string-string map type, defined in stdafx.h
	ar & p.m_mapUsed; // a string-bool map type, defined in stdafx.h
	ar & p.m_strXmlPath;
    }

  template<class Archive>
  void serialize(Archive &ar, mcleave_single &c, const unsigned int version)
  {
		ar & c.m_pNCleave;
		ar & c.m_pCCleave;
		ar & c.m_bN;
		ar & c.m_bC;
		ar & c.m_bNX;
		ar & c.m_bCX;
		ar & c.m_lType;
  }

  template<class Archive>
  void serialize(Archive &ar, mcleave &c, const unsigned int version)
  {
		ar & c.m_lType;
		ar & c.m_vCleaves;
		int s = (int)(c.m_itStart-c.m_vCleaves.begin());
		int e = (int)(c.m_itEnd-c.m_vCleaves.begin());
		ar & s;
		ar & e;
		c.m_itStart = c.m_vCleaves.begin()+s;
		c.m_itEnd = c.m_vCleaves.begin()+e;
	}

  template<class Archive>
  void serialize(Archive &ar, msequenceServer &p, const unsigned int version)
    {
  	ar & p.m_tColMax;		// Maximum size of a sequence collection
	ar & p.m_tStartAt;	// Ordinal number of the first sequence retreived
	ar & p.m_strPath;		// Full path name to the current FASTA file
	ar & p.m_strStatus;    // Status string for debugging
	ar & p.m_strFirst;
	ar & p.m_strTaxonPath;	// Path to the taxonomy translation file
	ar & p.m_strTaxon;		// Taxon to model
	ar & p.m_pCol;   // Sequence collection ring
	size_t initial_size = p.m_dstrFasta.size();
	size_t newsize = initial_size;
	ar & newsize;
	if (newsize != initial_size) { // we're writing
		p.m_dstrFasta.resize(newsize);
	}
	while (newsize--) {
		ar & p.m_dstrFasta[newsize];
	}
	ar & p.m_vstrFasta;
	ar & p.m_vstrDesc;
  }

  template<class Archive>
  void serialize(Archive &ar, merrors &p, const unsigned int version)
    {
	ar & p.m_bPpm;
	ar & p.m_bIsotope;
	ar & p.m_fPlus;
	ar & p.m_fMinus;
  }

template<class Archive>
  void serialize(Archive &ar,  mpyrostate &p, const unsigned int version)
    {
	ar & p.m_bPotential;
	ar & p.m_bPyro;
	ar & p.m_dModMass;
	ar & p.m_cRes;
  }

  template<class Archive>
  void serialize(Archive &ar,  msemistate  &p, const unsigned int version)
    {
  	  ar & p.m_bActive;
	  ar & p.m_lStart;
	  ar & p.m_lEnd;
	  ar & p.m_lStartI;
	  ar & p.m_lEndI;
	  ar & p.m_bStart;
	  ar & p.m_lLimit;
	  ar & p.m_lLastCleave;
  }

  template<class Archive>
  void serialize(Archive &ar, mprocess &p, const unsigned int version)
    {
      // omitted m_prcLog
      ar & p.m_notes;
      ar & p.m_hostname;
      ar & p.m_xmlPerformance;
      ar & p.m_xmlValues;
      ar & p.m_vSpectra;
	ar & p.m_mapSequences; // a map containing all of the protein sequences discovered, indexed by their m_tUid value
	ar & p.m_vseqBest; // a vector of msequences used in the model refinement process
	ar & p.m_vstrModifications; //a vector containing the strings defining fixed modifications for a protein (bpratt 9/1/2010)
	ar & p.m_tRefineModels; // total number of models generated by refinement
	ar & p.m_tRefineInput; // total number of sequences included in a refinement session
	ar & p.m_tRefinePartial; // the number of models discovered to have partial cleavage
	ar & p.m_tRefineUnanticipated; // the number of models discovered to have unanticpated cleavage
	ar & p.m_tRefineNterminal; // the number of models discovered to have modified N-terminii
	ar & p.m_tRefineCterminal; // the number of models discovered to have modified C-terminii
	ar & p.m_tRefinePam; // the number of models discovered to have point mutations
	ar & p.m_dRefineTime; // the time required to perform a refinement
	ar & p.m_tActive;	// total number of models remaining after each refinement step
	ar & p.m_bRefineCterm;  //true if processing 'refine, potential C-terminus modifications'. Set in mrefine::refine and 

	ar & p.m_viQuality; // contains the data quality scoring vector
	ar & p.m_bReversedOnly;
	ar & p.m_bSaps;
	ar & p.m_bAnnotation;

	ar & p.m_strLastMods;
	ar & p.m_iCurrentRound;
	ar & p.m_bPermute;
	ar & p.m_bPermuteHigh;
	ar & p.m_bCrcCheck;

	size_t original_n = (unsigned long)p.m_setRound.size();
	size_t n = original_n;
	ar & n;
	if (n != original_n) { // loading
		while (n--) {
			size_t r;
			ar & r;
			p.m_setRound.insert(r);
		}
	} else for (std::set<size_t>::iterator r = p.m_setRound.begin(); r!=p.m_setRound.end();r++ ) {
		size_t val = *r;
		ar & val;
	}
	ar & p.m_vstrSaps;
	ar & p.m_vstrMods;
	ar & p.m_mapAnnotation;

	ar & p.m_semiState; // maintains the state of the semi-enzymatic cleavage state machine
	ar & p.m_pyroState; // maintains the state of the pyrolidone carboxylic acid detection state machine
	ar & p.m_errValues;
	ar & p.m_dSearchTime; // total time elapsed during a protein modeling session process
	ar & p.m_lIonCount; // minimum sum of detected ions that are significant enough to store a sequence 
	ar & p.m_lThread; // thread number of this object
	ar & p.m_lThreads; // the total number of threads current active
	ar & p.m_lReversed; // the total number of peptides found where the reversed sequence was better than the forward sequence
	ar & p.m_dThreshold; // the current expectation value threshold
	ar & p.m_tContrasted; // the number of spectra subtracted using contrast angle redundancy detection
	ar & p.m_lStartMax; // set the maximum distance from N-terminus for a peptide
					  // normally set at an impossibly large value = 100,000,000
					  // for ragged N-terminii with potential modifications, set at a low but plausible value = 50
	ar & p.m_lCStartMax;
	// omitted ar & p.m_pSeq; // a character pointer, used for temporary sequence information
	ar & p.m_bUn; // if true, cleave at all residues. if false, use cleavage specification in data input.
	ar & p.m_bUseHomologManagement; // set to true to use homologue management 
	ar & p.m_tMinResidues; // the minimum peptide length that will be scored
	ar & p.m_tMissedCleaves; // the maximum number of cleavage sites that can be missed
	ar & p.m_tPeptideCount; // the total number of peptide sequences generated during a process
	ar & p.m_tPeptideScoredCount; // the total number of peptide sequences scored during a process
	ar & p.m_tProteinCount; // the total number of protein sequences considered during a process
	ar & p.m_tSpectra; // the total number of spectra being modeled
	ar & p.m_tSpectraTotal; // the total number of spectra in the input file
	ar & p.m_tValid; // the number of valid peptide models
	ar & p.m_tTotalResidues; // the number of residues read
	ar & p.m_tSeqSize; // current length of the m_pSeq character array
	ar & p.m_tUnique; // the number of unique peptides found in a result
	ar & p.m_strOutputPath; // the path name of the XML output file
	ar & p.m_Cleave; // the specification for a cleavage peptide bond
	ar & p.m_seqCurrent; // the msequence object that is currently being scored
	ar & p.m_svrSequences; // the msequenceServer object that provides msequences to msequenceCollection
	ar & p.m_specCondition; // the mspectrumcondition object that cleans up and normalized
										// spectra for further processing
	ar & p.m_pScore; // the object that is used to score sequences and spectra
	// omitted mrefine* m_pRefine; // the object that is used to refine models
    }

  template<class Archive>
  void serialize(Archive &ar, maa &m, const unsigned int version)
    {
	ar & m.m_lPos; // the sequence position of the residue (N-terminal = 0)
	ar & m.m_dMod; // mass of the modification
	ar & m.m_cRes; // single letter abbreviation for the amino acid
	ar & m.m_cMut; // single letter abbreviation for a discovered point mutation
	ar & m.m_strId; // character string representing an external accession number for a mutation/modification
        ar & m.m_dPrompt; // prompt loss from modification mass
    }
  template<class Archive>
  void serialize(Archive &ar, mspectrum &s, const unsigned int version)
    {

	ar & s.m_tId; // an identification number
	ar & s.m_tCurrentSequence; // an identifier for the current sequence (used in mprocess)
	ar & s.m_fScore; // the convolution score
	ar & s.m_fHyper; // the hyper score
	ar & s.m_fScoreNext; // next best convolution score
	ar & s.m_fHyperNext; // next best hyper score
	ar & s.m_dExpect; // the expectation value
	ar & s.m_dProteinExpect; // the expectation value for the associated protein
	ar & s.m_dMH; // the parent ion mass + a proton
	ar & s.m_fI; // the parent ion intensity (if available)
	ar & s.m_fZ; // the parent ion charge
	ar & s.m_bRepeat; // a flag indicating that a better match for an individual peptide has already been found
	ar & s.m_bActive; // a flag indicating that a spectrum is available for scoring
	ar & s.m_vMI; // a vector containing the m/z - intensity information for fragment ions 
	ar & s.m_vMINeutral; // a vector containing the m/z - intensity information for fragment ions 
	ar & s.m_vseqBest; // a vector containing the highest scoring msequence objects
	ar & s.m_vdStats;
	ar & s.m_strDescription;
 	ar & s.m_strRt;

	ar & s.m_hHyper; // the histogram of hyper scores
	ar & s.m_hConvolute; // the histogram of convolution scores
	ar & s.m_chBCount; // the histogram of b-ion counts
	ar & s.m_chYCount; // the histogram of y-ion counts
    }
  template<class Archive>
  void serialize(Archive &ar, mi &m, const unsigned int version)
    {
  	ar & m.m_fM; // the m/z value
        ar & m.m_fI; // the intensity value
    }
  template<class Archive>
  void serialize(Archive &ar, mdomain &d, const unsigned int version)
    {
	ar & d.m_lS; // the start position of the peptide in the protein sequence (N-terminus = 0)
	ar & d.m_lE; // the end position of the peptide in the protein sequence
	ar & d.m_lMissedCleaves; // missed cleavages
	ar & d.m_fScore; // the convolution score for the peptide
        ar & d.m_fHyper; // the hyper score for the peptide 
	ar & d.m_dMH; // the mass of the peptide + a proton
	// double m_dDelta replaces float m_fDelta, starting with version 2006.02.01
	// because of an issue with the accuracy of this value
	ar & d.m_dDelta; // the mass difference between the mass of the peptide and the measured mass
	ar & d.m_bUn;
        ar & d.m_mapCount; // a map of the number of ions detected for each ion type
	ar & d.m_mapScore; // a map of the convolution scores for each ion type
	ar & d.m_vAa; // vector of modified amino acids
    }


  template<class Archive>
  void serialize(Archive &ar, mspectrumcondition &s, const unsigned int version)
    {
	ar & s.m_bCondition; // enables the use of all conditioning methods
	ar & s.m_bUseChargeSuppression; // enables the rejection of highly charge parent ions
	ar & s.m_bUseDynamicRange; // enables using the dynamic range
	ar & s.m_bUseLowestMass; // enables the removal of very low m/z peaks from a spectrum
	ar & s.m_bUseMaxPeaks; // enables removing low intensity peaks
	ar & s.m_bUseMinMass; // sets the minimum parent ion mass allowed
	ar & s.m_bUseMaxMass; // sets the maximum parent ion mass allowed
	ar & s.m_bUseMinSize; // enables using the minimum number of peaks to exclude spectra
	ar & s.m_bUseNoiseSuppression; // enables the detection of purely noise spectra
	ar & s.m_bUseParent; // enables the exclusion of spectra by the parent ion mass
	ar & s.m_bUseNeutralLoss;
	ar & s.m_bUsePhosphoDetection;
	ar & s.m_tMaxPeaks; // the maximum number of peaks in a spectrum
	ar & s.m_fDynamicRange; // the normalized intensity of the most intense peak in a spectrum
	ar & s.m_fLowestMass; // the lowest m/z in a spectrum
	ar & s.m_lMinSize; // the minimum number of peaks in a spectrum
	ar & s.m_fMinMass; // the minimum parent ion mass in a spectrum
	ar & s.m_fMaxMass; // the maximum parent ion mass in a spectrum
	ar & s.m_fParentLower; // the low end of the mass window for excluding parent ion neutral losses
	ar & s.m_fParentUpper; // the high end of the mass window for excluding parent ion neutral losses (just passed the last C13 isotope peak)
	ar & s.m_fMaxZ; // the maximum parent ion charge allowed
	ar & s.m_fNeutralLoss;
	ar & s.m_fNeutralLossWidth;
	ar & s.m_fFactor;
	}



}
}

// save process to a compressed, base64'd string (caller must delete[])
// result is formatted as 
// <length of base64 string><tab><length of data after decompressing><tab><base64'ed data>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream_buffer.hpp> 
#include <boost/iostreams/stream.hpp> 
#include <boost/iostreams/device/back_inserter.hpp> 
#include <boost/serialization/vector.hpp> 
#include <vector> 
#include <iostream> 
namespace io = boost::iostreams;

// result is formatted as 
// <length of base64 string><tab><length of data after decompressing><tab><base64'ed data>
// or
// <length of compressed data><tab><length of data after decompressing><tab><compressed data>
char *serialize_process(mprocess &p, int *resultlen, bool bAsText){
  std::vector<char> outputbuffer;
  io::stream<io::back_insert_device< std::vector<char> > > oss(outputbuffer); 
  oarchive oa(oss,ios::binary);
  oa << p;
  oss.flush();
  unsigned long olen = (unsigned long)outputbuffer.size();
  const char *obuf = &outputbuffer[0];
  unsigned long zbuflen;
  Bytef *zbuf=new Bytef [zbuflen=compressBound(olen)+1];
  compress(zbuf,&zbuflen,(const Bytef *)obuf,olen);
  size_t bbuflen;
  char *result = new char[100+(bbuflen=pwiz::util::Base64::binaryToTextSize(zbuflen))];
  sprintf(result,"%u\t%u\t",bAsText?bbuflen:zbuflen,olen);
  int rlen=strlen(result);
  char *bbuf = result+rlen;
  if (bAsText) { // base64 it
    size_t blen=pwiz::util::Base64::binaryToText((char *)zbuf,(int)zbuflen, bbuf);
    bbuf[blen] = 0;
    *resultlen = rlen+blen;
  } else {
	memmove(bbuf,zbuf,zbuflen);
    *resultlen = rlen+zbuflen;
  }
  delete[] zbuf;
  return result;
}

void deserialize_process(int filehandle, bool bAsText, mprocess &p,bool verbose){
	// read the raw data length
	const int MAXHDRSIZE=512;
	char buf[MAXHDRSIZE];
	int nread=read(filehandle,buf,MAXHDRSIZE);
	int rawlen = atoi(buf)+MAXHDRSIZE; // allow some extra 
	char *databuf = new char[rawlen+1];
	memmove(databuf,buf,nread);
	int nreadnext = read(filehandle,databuf+nread,rawlen-nread);
	deserialize_process(databuf, bAsText, p, verbose);
	delete[] databuf;
}

void deserialize_process(const char *data, bool bAsText, mprocess &p,bool verbose){
	size_t raw_len = (size_t)atol(data);
	size_t uncompressed_len = (size_t)atol(strchr(data,'\t')+1);
	const char *encodedData = strchr(strchr(data,'\t')+1,'\t')+1;
	Bytef *processdata = new Bytef[uncompressed_len+1];
	uLongf processdataLen = uncompressed_len+1;
	if (bAsText) {
		char *compressedData = new char[raw_len];
		uLong compressedLen = (uLong)pwiz::util::Base64::textToBinary(encodedData, raw_len, compressedData);
	uncompress(processdata,&processdataLen,(Bytef *)compressedData,compressedLen);
	delete[] compressedData;
	} else {
		uncompress(processdata,&processdataLen,(Bytef *)encodedData,raw_len);
	}
	io::basic_array_source<char> source((char *)processdata,processdataLen);
	io::stream<io::basic_array_source <char> > iss(source);
	iarchive ia(iss,ios::binary);
	ia >> p;

	delete[] processdata;
	if (verbose) {
		std::cerr << mapreducehelper::timestamp() << "loaded process with spectra count = " << (int)p.m_vSpectra.size() << "\n";
	}
}

#include "mscore_hrk.h"
#include "mscore_c.h"
#include "mscore_k.h"
#include "mscore_tandem.h"
BOOST_CLASS_EXPORT(mscore_k); // needed for serialization
BOOST_CLASS_EXPORT(mscore_hrk); // needed for serialization
BOOST_CLASS_EXPORT(mscore_c); // needed for serialization
BOOST_CLASS_EXPORT(mscore_tandem); // needed for serialization


  // our serialized objects are onetime use, don't bother with version info
#define NO_VERSIONING(classname) BOOST_CLASS_IMPLEMENTATION(classname, boost::serialization::object_serializable)
#define NO_PTR_TRACKING(classname) BOOST_CLASS_TRACKING(classname, boost::serialization::track_never)
#define MINIMAL_OVERHEAD_SERIALIZATION(classname) NO_VERSIONING(classname);NO_PTR_TRACKING(classname)
//
/* these need to be fully serialized, apparently because pointers are involved
  NO_VERSIONING(msequenceCollection);
  NO_VERSIONING(msequenceServer);
  NO_VERSIONING(mscore_c);
  NO_VERSIONING(mscore_hrk);
  NO_VERSIONING(mscore_k);
  NO_VERSIONING(mscore_tandem);
  NO_VERSIONING(mscore);
*/
  MINIMAL_OVERHEAD_SERIALIZATION(XmlParameter);
  MINIMAL_OVERHEAD_SERIALIZATION(mcleave_single);
  MINIMAL_OVERHEAD_SERIALIZATION(mcleave);
  MINIMAL_OVERHEAD_SERIALIZATION(merrors);
  MINIMAL_OVERHEAD_SERIALIZATION(mpyrostate);
  MINIMAL_OVERHEAD_SERIALIZATION(msemistate);
  MINIMAL_OVERHEAD_SERIALIZATION(maa);
  MINIMAL_OVERHEAD_SERIALIZATION(mspectrum);
  MINIMAL_OVERHEAD_SERIALIZATION(mi);
  MINIMAL_OVERHEAD_SERIALIZATION(mspectrumcondition);
  MINIMAL_OVERHEAD_SERIALIZATION(mdomain);
  MINIMAL_OVERHEAD_SERIALIZATION(PermuteState);
  MINIMAL_OVERHEAD_SERIALIZATION(mspectrumindex);
  MINIMAL_OVERHEAD_SERIALIZATION(mspec);
  MINIMAL_OVERHEAD_SERIALIZATION(MIType);
  MINIMAL_OVERHEAD_SERIALIZATION(mspectrumdetails);
  MINIMAL_OVERHEAD_SERIALIZATION(masscalc); 
  MINIMAL_OVERHEAD_SERIALIZATION(masscalc::massPair); 
  MINIMAL_OVERHEAD_SERIALIZATION(mhistogram); 
  MINIMAL_OVERHEAD_SERIALIZATION(count_mhistogram); 
  MINIMAL_OVERHEAD_SERIALIZATION(mmotifres);
  MINIMAL_OVERHEAD_SERIALIZATION(mmotif);
  MINIMAL_OVERHEAD_SERIALIZATION(msequence);
  MINIMAL_OVERHEAD_SERIALIZATION(msequtilities);

#endif // HAVE_MULTINODE_TANDEM
