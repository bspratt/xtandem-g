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

// File version: 2004-01-07
// File version: 2005-11-15

/*
 * loadspectrum.cpp contains the override methods necessary for compatibility with
 * various mass spectrum file formats. See loadspectrum.h for additional details.
 */

#include "stdafx.h"
#include <algorithm>
#include "msequence.h"
#include "mspectrum.h"
#include "loadmspectrum.h"
#include "base64.h"

loadgaml::loadgaml( vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m)
	: handler(_vS, _sC, _m)
{
}

loadgaml::~loadgaml(void)
{
}

bool loadgaml::get()
{
  handler.parse();
  return true;
}

bool loadgaml::open(string &_s)
{
	/*
	* copy the file pathname into m_strPath
	*/
	m_strPath = _s;
	/*
	* try to open the file and bail out if it isn't available
	*/
	m_ifIn.open(m_strPath.c_str());
	if(m_ifIn.fail())	{
		return false;
	}

  /*
  * test for the file extension .bioml
  */
  string strTest = m_strPath;
  int (*pf)(int) = tolower; 
  transform(strTest.begin(), strTest.end(), strTest.begin(), pf); 
  if(strTest.find(".bioml") != strTest.npos)	{
	  m_ifIn.close();
	  handler.setFileName( m_strPath.c_str() );
	  return true;
  }
  /*
   * check the open file to see if it is a BIOML file
   * Verification tres tres sommaire
   * scoop up the first 32 lines and then test for
   * <?xml followed by <bioml
   */
  long a = 0;
  strTest =	" ";
  string strLine;
  getline(m_ifIn,strLine);
  strTest += strLine;
  while(m_ifIn.good() && !m_ifIn.eof() && a	< 32)	{
	  getline(m_ifIn,strLine);
	  strTest += strLine;
	  a++;
  }
  m_ifIn.close();
  size_t tXml =	strTest.find("<?xml");
  size_t tMz = 0;
  if(tXml != strTest.npos)	{
	  tMz =	strTest.find("xmlns:GAML=",tXml);
	  if(tMz ==	strTest.npos)	{
			return false;
	  }
  }
  else	{
	  return false;
  }

	handler.setFileName( m_strPath.c_str() );
	return true;
}
//loadgaml::open_force
// forces the opening of an GAML file
bool loadgaml::open_force(string &_s)
{
	/*
	* copy the file pathname into m_strPath
	*/
	m_strPath = _s;
	/*
	* try to open the file and bail out if it isn't available
	*/
	m_ifIn.open(m_strPath.c_str());
	if(m_ifIn.fail())	{
		return false;
	}

	m_ifIn.close();
	m_ifIn.clear();

	handler.setFileName( m_strPath.c_str() );

	return true;
}

#ifdef XMLCLASS
/*
loadmzdata est ajouter par Patrick Lacasse en janvier 2005
 */
loadmzdata::loadmzdata( vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m)
	: handler(_vS, _sC, _m)
{
}

loadmzdata::~loadmzdata(void)
{
}

/*
 * 
 * loadmzdata::get
 * loadmzdata::open doit avoir ete lancee avant.
 * parser est SAX2XMLparser* instantie par le constraucteur
 * get met tous les spectres du fichier dans _vSpec
 *
 */
bool loadmzdata::get()
{
  handler.parse();
  return true;
}


// loadmzdata::open
// La machine xml est ouverte ici en meme temps que le fichier d'entree.
// Elle sera refermee lorsque loadmzdata::get atteindra la fin de m_ifIn
bool loadmzdata::open(string &_s)
{
  //Il faut trouver une fonction propre pour tester qu'un fichier est un mzXML
  m_tId = 1;
  /*
   * copy the file pathname into m_strPath
   */
  m_strPath = _s;
  /*
   * try to open the file and bail out if it isn't available
   */
  m_ifIn.open(m_strPath.c_str());
  if(m_ifIn.fail())	{
    return false;
  }
  /*
  * test for the file extension .mzdata
  */
  string strTest = m_strPath;
  int (*pf)(int) = tolower; 
  transform(strTest.begin(), strTest.end(), strTest.begin(), pf); 
  if(strTest.find(".mzdata") != strTest.npos)	{
	  m_ifIn.close();
	  handler.setFileName( m_strPath.c_str() );
	  return true;
  }
  /*
   * check the open file to see if it is an mzData file
   * Verification tres tres sommaire
   * scoop up the first 32 lines and then test for
   * <?xml followed by <mzData
   */
  long a = 0;
  strTest =	" ";
  getline(m_ifIn,strLine);
  strTest += strLine;
  while(m_ifIn.good() && !m_ifIn.eof() && a	< 32)	{
	  getline(m_ifIn,strLine);
	  strTest += strLine;
	  a++;
  }
  m_ifIn.close();
  size_t tXml =	strTest.find("<?xml");
  size_t tMz = 0;
  if(tXml != strTest.npos)	{
	  tMz =	strTest.find("<mzData",tXml);
	  if(tMz ==	strTest.npos)	{
			return false;
	  }
  }
  else	{
	  return false;
  }
  handler.setFileName( m_strPath.c_str() );

  return true;
}
// loadmzdata::open_force
// forces the opening of an mzDATA file
bool loadmzdata::open_force(string &_s)
{
  //Il faut trouver une fonction propre pour tester qu'un fichier est un mzXML
  m_tId = 1;
  /*
   * copy the file pathname into m_strPath
   */
  m_strPath = _s;
  /*
   * try to open the file and bail out if it isn't available
   */
  m_ifIn.open(m_strPath.c_str());
  if(m_ifIn.fail())	{
    return false;
  }
  m_ifIn.close();
  handler.setFileName( m_strPath.c_str() );
  return true;
}


/*
loadmzxml est ajouter par Patrick Lacasse en decembre 2004
 */
loadmzxml::loadmzxml( vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m)
	: handler(_vS, _sC, _m)
{
}

loadmzxml::~loadmzxml(void)
{
}

/*
 * 
 * loadmzxml::get
 * loadmzxml::open doit avoir ete lancee avant.
 * parser est SAX2XMLparser* instantie par le constraucteur
 * get met tous les spectres du fichier dans _vSpec
 *
 */
bool loadmzxml::get()
{  
  return handler.parse();
}


// loadmzxml::open
// La machine xml est ouverte ici en meme temps que le fichier d'entree.
// Elle sera refermee lorsque loadmzxml::get atteindra la fin de m_ifIn
bool loadmzxml::open(string &_s)
{
  //Il faut trouver une fonction propre pour tester qu'un fichier est un mzXML
  m_tId = 1;
  /*
   * copy the file pathname into m_strPath
   */
  m_strPath = _s;
	  
  /*
   * try to open the file and bail out if it isn't available
   */
  m_ifIn.open(m_strPath.c_str());
  if(m_ifIn.fail())	{
    return false;
  }
  /*
  * test for the file extension .mzxml
  */
  string strTest = m_strPath;
  int (*pf)(int) = tolower; 
  transform(strTest.begin(), strTest.end(), strTest.begin(), pf); 
  if(strTest.find(".mzxml") != strTest.npos)	{
	  m_ifIn.close();
	  handler.setFileName( m_strPath.c_str() );
	  return true;
  }
  /*
   * check the open file to see if it is an mzXML file
   * Verification tres tres sommaire
    * scoop up the first 32 lines and then test for
   * <?xml followed by <mzXML or <msRun
   */
  long a = 0;
  strTest = " ";
  getline(m_ifIn,strLine);
  strTest += strLine;
  while(m_ifIn.good() && !m_ifIn.eof() && a	< 32)	{
	  getline(m_ifIn,strLine);
	  strTest += strLine;
	  a++;
  }
  m_ifIn.close();
  size_t tXml = strTest.find("<?xml");
  size_t tMz = 0;
  if(tXml != strTest.npos)	{
	  // Some version of ReadW.exe create a	file with msRun	as the primary tag,
	  // rather	than mzXML,	so each	should be taken	into account. -S.Wiley
	  tMz = strTest.find("<mzXML",tXml);
	  if(tMz == strTest.npos)	{
		  tMz = strTest.find("<msRun",tXml);
		  if(tMz == strTest.npos)	{
			  return false;
		  }
	  }
  }
  else	{
	  return false;
  }
  handler.setFileName( m_strPath.c_str() );
  return true;
}
// loadmzxml::open_force
// forces the opening of an mzXML file
bool loadmzxml::open_force(string &_s)
{
  //Il faut trouver une fonction propre pour tester qu'un fichier est un mzXML
  m_tId = 1;
  /*
   * copy the file pathname into m_strPath
   */
  m_strPath = _s;
  /*
   * try to open the file and bail out if it isn't available
   */
  m_ifIn.open(m_strPath.c_str());
  if(m_ifIn.fail())	{
    return false;
  }
  m_ifIn.close();
  handler.setFileName( m_strPath.c_str() );
  return true;
}

#endif //ifdef XMLCLASS

loadmatrix::loadmatrix(void)
{
}

loadmatrix::~loadmatrix(void)
{
}

bool loadmatrix::get(mspectrum &_m)
{
	char *pLine = new char[m_tSize];
	char *pValue;
	bool bFirst = true;
	long lId = 0;
	mspectrum specCurrent;
	specCurrent.m_strDescription.erase(specCurrent.m_strDescription.begin(),specCurrent.m_strDescription.end());
/*
 * find the next spectrum in the file - the stream should still be connected
 */
	while(!m_ifIn.eof() && m_ifIn.good())	{
		m_ifIn.getline(pLine,m_tSize-1,m_cEol);
		if(strstr(pLine,"BEGIN IONS") != NULL)
			break;
	}
	mi miCurrent;
	bool bNext = true;
/*
 * create a temporary mspectrum object
 */
	specCurrent.clear_intensity_values();
	specCurrent.m_fZ = 2.0;
/*
 * read the parent ion mass and charge and find all of the fragment m/z values and intensities
 */
	while(m_ifIn.good() && !m_ifIn.eof())	{
		m_ifIn.getline(pLine,m_tSize-1,m_cEol);
		if(strstr(pLine,"PEPMASS=") != NULL)	{
			pValue = strchr(pLine,'=');
			pValue++;
			specCurrent.m_dMH = atof(pValue);
		}
		else if(strstr(pLine,"#") == pLine)	{
			specCurrent.m_strDescription += pLine+1;
			specCurrent.m_strDescription += " ";
		}
		else if(strstr(pLine,"TITLE=") == pLine)	{
			specCurrent.m_strDescription += strchr(pLine,'=')+1;
			specCurrent.m_strDescription += " ";
		}
		else if(strstr(pLine,"CHARGE=") != NULL)	{
			pValue = strchr(pLine,'=');
			pValue++;
			specCurrent.m_fZ = (float)atof(pValue);
		}
		else if(atof(pLine) > 0.0)	{
			miCurrent.m_fM = (float)atof(pLine);
			pValue = pLine;
			while(*pValue != '\0' && isspace(*pValue))	{
				pValue++;
			}
			while(*pValue != '\0' && !isspace(*pValue))	{
				pValue++;
			}
			miCurrent.m_fI = (float)atof(pValue);
			specCurrent.m_vMI.push_back(miCurrent);
		}
		else if(strstr(pLine,"END IONS") != NULL)	{
			break;
		}
	}
	delete pLine;
	double dProton = 1.007276;
/*
 * adjust the parent ion mass to be M+H (mspectrum expects an M+H) 
 */

	specCurrent.m_dMH = (specCurrent.m_dMH - dProton)*specCurrent.m_fZ + dProton;
	specCurrent.m_tId = m_tId;
/*
 * copy the temporary mspectrum into the input reference 
 */
	_m = specCurrent;
	m_tId++;
/*
 * if the file is finished, return false, otherwise return true 
 */
	if(m_ifIn.eof() || !m_ifIn.good())	{
		m_ifIn.close();
		return false;
	}
	return true;
}

bool loadmatrix::open(string &_s)
{
	m_tId = 1;
/*
 * copy the file pathname into m_strPath
 */
	m_strPath = _s;
/*
 * try to open the file and bail out if it isn't available
 */
	m_ifIn.open(m_strPath.c_str());
	if(m_ifIn.fail())	{
		return false;
	}
/*
 * check the open file to see if it is a matrix science file
 */
	char *pLine = new char[m_tSize];
	m_ifIn.getline(pLine,256,m_cEol);
	pLine[255] = '\0';
	if(strlen(pLine) == 255)	{
		m_cEol = 0x0D;
	}
	m_ifIn.close();
	m_ifIn.clear();
	m_ifIn.open(m_strPath.c_str());
	m_ifIn.getline(pLine,m_tSize-1,m_cEol);
	pLine[m_tSize-1] = '\0';
	bool bOk = false;
	long lCount = 0;
	while(!bOk && !m_ifIn.eof() && lCount < 4096)	{
		if(strstr(pLine,"BEGIN IONS") == pLine)	{
			bOk = true;
		}
		lCount++;
		m_ifIn.getline(pLine,m_tSize-1,m_cEol);
		pLine[m_tSize-1] = '\0';
	}
	m_ifIn.close();
/*
 * return false if there were no matrix science spectrum tags
 */
	if(!bOk)	{
		delete pLine;
		return false;
	}
/*
 * if it is a matrix science file, clear the stream and reopen it for calls from get
 */
	m_ifIn.clear();
	m_ifIn.open(m_strPath.c_str());
	delete pLine;
	return true;
}
// loadmatrix::open_force
// forces the opening of a Matrix Science Generic format file
bool loadmatrix::open_force(string &_s)
{
	m_tId = 1;
/*
 * copy the file pathname into m_strPath
 */
	m_strPath = _s;
/*
 * try to open the file and bail out if it isn't available
 */
	m_ifIn.open(m_strPath.c_str());
	if(m_ifIn.fail())	{
		return false;
	}
	m_ifIn.close();
/*
 * return false if there were no matrix science spectrum tags
 */
	m_ifIn.clear();
	m_ifIn.open(m_strPath.c_str());
	return true;
}

loadpkl::loadpkl(void)
{
}

loadpkl::~loadpkl(void)
{
}


bool loadpkl::get(mspectrum &_m)
{
	char *pLine = new char[m_tSize];
	char *pValue;
	bool bFirst = true;
	long lId = 0;
	mi miCurrent;
	bool bNext = true;
/*
 * create a temporary mspectrum object
 */
	mspectrum specCurrent;
	specCurrent.m_strDescription = "no description";
	double dProton = 1.007276;
	specCurrent.m_fZ = 2.0;
	while(bNext && m_ifIn.good() && !m_ifIn.eof())	{
		m_ifIn.getline(pLine,m_tSize-1,m_cEol);
/*
 * find the next spectrum in the file - the stream should still be connected
 */
		if(atof(pLine) != 0.0)	{
			if(bFirst)	{
				specCurrent.clear_intensity_values();
/*
 * load the parent m/z value 
 */
				specCurrent.m_dMH = atof(pLine);
				pValue = pLine;
				while(*pValue != '\0' && isspace(*pValue))	{
					pValue++;
				}
				while(*pValue != '\0' && !isspace(*pValue))	{
					pValue++;
				}
/*
 * find and skip the parent ion intensity 
 */
				while(*pValue != '\0' && isspace(*pValue))	{
					pValue++;
				}
				while(*pValue != '\0' && !isspace(*pValue))	{
					pValue++;
				}
/*
 * load the parent charge value 
 */
				if(*pValue != '\0')	{
					specCurrent.m_fZ = (float)atof(pValue);
				}
/*
 * convert the parent m/z value into M+H, as required by mspectrum definition 
 */
				specCurrent.m_dMH = (specCurrent.m_dMH - dProton)*specCurrent.m_fZ + dProton;
				bFirst = false;
			}
			else	{
/*
 * load fragment ion m/z values and intensities 
 */
				miCurrent.m_fM = (float)atof(pLine);
				pValue = pLine;
				while(*pValue != '\0' && isspace(*pValue))	{
					pValue++;
				}
				while(*pValue != '\0' && !isspace(*pValue))	{
					pValue++;
				}
				miCurrent.m_fI = (float)atof(pValue);
				while(*pValue != '\0' && isspace(*pValue))	{
					pValue++;
				}
				while(*pValue != '\0' && !isspace(*pValue))	{
					pValue++;
				}
				if(pValue != '\0' && strlen(pValue) > 2)	{
					specCurrent.m_strDescription = pValue;
				}
				specCurrent.m_vMI.push_back(miCurrent);
			}
		}
		else	{
			if(specCurrent.m_vMI.size() > 0)	{
				bNext = false;
			}
			bFirst = true;
		}
	}
	delete pLine;
	specCurrent.m_tId = m_tId;
/*
 * copy the temporary mspectrum into the input reference 
 */
	_m = specCurrent;
	m_tId++;
/*
 * if the file is finished, return false, otherwise return true 
 */
	if(m_ifIn.eof())	{
		m_ifIn.close();
		return false;
	}
	return true;
}

bool loadpkl::open(string &_s)
{
	m_tId = 1;
/*
 * copy the file pathname into m_strPath
 */
	m_strPath = _s;
/*
 * try to open the file and bail out if it isn't available
 */
	m_ifIn.open(m_strPath.c_str());
	if(m_ifIn.fail())	{
		return false;
	}
/*
 * check the open file to test for the value of the end-of-line character
 */
	long lCount = 0;
	char *pLine = new char[m_tSize];
	m_ifIn.getline(pLine,256,m_cEol);
	pLine[255] = '\0';
	if(strlen(pLine) == 255)	{
		m_cEol = 0x0D;
	}
	m_ifIn.close();
	m_ifIn.clear();
/*
 * close, clear and reopen the stream to test to see if it is a PKL file
 * the test is to see if the first non-blank line has three ASCII numbers separated
 * by whitespace
 */
	m_ifIn.open(m_strPath.c_str());
	m_ifIn.getline(pLine,m_tSize-1,m_cEol);
	while(!m_ifIn.eof() && atof(pLine) == 0.0 && lCount < 4096)	{
		m_ifIn.getline(pLine,m_tSize-1,m_cEol);
		lCount++;
	}
	if(m_ifIn.eof())	{
		m_ifIn.close();
		delete pLine;
		return false;
	}
	char *pValue = pLine;
	while(*pValue != '\0' && isspace(*pValue))	{
		pValue++;
	}
	while(*pValue != '\0' && !isspace(*pValue))	{
		pValue++;
	}
	if(atof(pValue) == 0.0)	{
		m_ifIn.close();
		delete pLine;
		return false;
	}
	while(*pValue != '\0' && isspace(*pValue))	{
		pValue++;
	}
	while(*pValue != '\0' && !isspace(*pValue))	{
		pValue++;
	}
	if(atof(pValue) == 0.0)	{
		m_ifIn.close();
		delete pLine;
		return false;
	}
	m_ifIn.close();
	m_ifIn.clear();
/*
 * if the file format is ok, clear, close and reopen the stream, returning true
 */
	m_ifIn.open(m_strPath.c_str());
	delete pLine;
	return true;
}
// loadpkl::open_force
// forces the opening of a PKL format file
bool loadpkl::open_force(string &_s)
{
	m_tId = 1;
/*
 * copy the file pathname into m_strPath
 */
	m_strPath = _s;
/*
 * try to open the file and bail out if it isn't available
 */
	m_ifIn.open(m_strPath.c_str());
	if(m_ifIn.fail())	{
		return false;
	}
/*
 * check the open file to test for the value of the end-of-line character
 */
	char *pLine = new char[m_tSize];
	m_ifIn.getline(pLine,256,m_cEol);
	pLine[255] = '\0';
	if(strlen(pLine) == 255)	{
		m_cEol = 0x0D;
	}
	m_ifIn.close();
	m_ifIn.clear();
/*
 * if the file format is ok, clear, close and reopen the stream, returning true
 */
	m_ifIn.open(m_strPath.c_str());
	delete pLine;
	return true;
}

loaddta::loaddta(void)
{
}

loaddta::~loaddta(void)
{
}

bool loaddta::get(mspectrum &_m)
{
	char *pLine = new char[m_tSize];
	char *pValue;
	bool bFirst = true;
	long lId = 0;
	mi miCurrent;
	bool bNext = true;
/*
 * create a temporary mspectrum object 
 */
	mspectrum specCurrent;
	specCurrent.m_strDescription = "no description";
	while(bNext && m_ifIn.good() && !m_ifIn.eof())	{
		m_ifIn.getline(pLine,m_tSize-1,m_cEol);
/*
 * find the next spectrum  
 */
		if(atof(pLine) != 0.0)	{
			if(bFirst)	{
				specCurrent.clear_intensity_values();
/*
 *  load the parent M+H value
 */
				specCurrent.m_dMH = atof(pLine);
				pValue = pLine;
				while(*pValue != '\0' && isspace(*pValue))	{
					pValue++;
				}
				while(*pValue != '\0' && !isspace(*pValue))	{
					pValue++;
				}
/*
 *  load the parent charge value
 */
				if(*pValue != '\0')	{
					specCurrent.m_fZ = (float)atof(pValue);
				}
				while(*pValue != '\0' && isspace(*pValue))	{
					pValue++;
				}
				while(*pValue != '\0' && !isspace(*pValue))	{
					pValue++;
				}
				if(pValue != '\0' && strlen(pValue) > 2)	{
					specCurrent.m_strDescription = pValue;
				}
				bFirst = false;
			}
			else	{
/*
 * load the m/z and intensity values  
 */
				miCurrent.m_fM = (float)atof(pLine);
				pValue = pLine;
				while(*pValue != '\0' && isspace(*pValue))	{
					pValue++;
				}
				while(*pValue != '\0' && !isspace(*pValue))	{
					pValue++;
				}
				miCurrent.m_fI = (float)atof(pValue);
				specCurrent.m_vMI.push_back(miCurrent);
			}
		}
		else	{
			if(specCurrent.m_vMI.size() > 0)	{
				bNext = false;
			}
			bFirst = true;
		}
	}
	delete pLine;
	specCurrent.m_tId = m_tId;
/*
 * copy the temporary mspectrum into the input reference
 */
	_m = specCurrent;
	m_tId++;
/*
 * if the file is finished, return false, otherwise return true 
 */
	if(m_ifIn.eof())	{
		m_ifIn.close();
		return false;
	}
	return true;
}

bool loaddta::open(string &_s)
{
	m_tId = 1;
/*
 * copy the file pathname into m_strPath
 */
	m_strPath = _s;
/*
 * try to open the file and bail out if it isn't available
 */
	m_ifIn.open(m_strPath.c_str());
	if(m_ifIn.fail())	{
		return false;
	}
/*
 * check the open file to test for the value of the end-of-line character
 */
	char *pLine = new char[m_tSize];
	m_ifIn.getline(pLine,256,m_cEol);
	pLine[255] = '\0';
	if(strlen(pLine) == 255)	{
		m_cEol = 0x0D;
	}
	m_ifIn.close();
	m_ifIn.clear();
/*
 * close, clear and reopen the stream to test to see if it is a DTA file
 * the test is to see if the first non-blank line has two ASCII numbers separated
 * by whitespace
 * note: PKL files also pass this test, so test for PKL files before testing for DTA files
 */
	m_ifIn.open(m_strPath.c_str());
	m_ifIn.getline(pLine,m_tSize-1,m_cEol);
	while(!m_ifIn.eof() && atof(pLine) == 0.0)	{
		m_ifIn.getline(pLine,m_tSize-1,m_cEol);
	}
	if(m_ifIn.eof())	{
		m_ifIn.close();
		delete pLine;
		return false;
	}
	char *pValue = pLine;
	while(*pValue != '\0' && isspace(*pValue))	{
		pValue++;
	}
	while(*pValue != '\0' && !isspace(*pValue))	{
		pValue++;
	}
	double dZ = atof(pValue);
	if(dZ == 0.0 || ((double)(int)dZ != dZ))	{
		m_ifIn.close();
		delete pLine;
		return false;
	}
	m_ifIn.close();
/*
 * if the file format is ok, clear, close and reopen the stream, returning true
 */
	m_ifIn.clear();
	m_ifIn.open(m_strPath.c_str());
	delete pLine;
	return true;
}
// loaddta::open_force
// forces the opening of a DTA format file
bool loaddta::open_force(string &_s)
{
	m_tId = 1;
/*
 * copy the file pathname into m_strPath
 */
	m_strPath = _s;
/*
 * try to open the file and bail out if it isn't available
 */
	m_ifIn.open(m_strPath.c_str());
	if(m_ifIn.fail())	{
		return false;
	}
/*
 * check the open file to test for the value of the end-of-line character
 */
	char *pLine = new char[m_tSize];
	m_ifIn.getline(pLine,256,m_cEol);
	pLine[255] = '\0';
	if(strlen(pLine) == 255)	{
		m_cEol = 0x0D;
	}
	m_ifIn.close();
	m_ifIn.clear();
	m_ifIn.open(m_strPath.c_str());
	delete pLine;
	return true;
}

loadcmn::loadcmn(void)
{
}

loadcmn::~loadcmn(void)
{
}

bool loadcmn::get(mspectrum &_m)
{
	if(m_pFile == NULL || feof(m_pFile))	{
		return false;
	}
	char *pLine = new char[256];
	long lId = 0;
	mi miCurrent;
	bool bNext = true;
/*
 * create a temporary mspectrum object 
 */
	mspectrum specCurrent;
	specCurrent.m_strDescription = "no description";
	unsigned short sValue = 0;
	unsigned char cValue = 0;
	unsigned int iValue = 0;
	float fValue = 0.0;
	double dValue = 0.0;

	size_t nread = fread((void *)&iValue,sizeof(unsigned int),1,m_pFile);
	m_tId = iValue;
	nread = fread((void *)&dValue,sizeof(double),1,m_pFile);
	specCurrent.m_dMH = dValue;
	nread = fread((void *)&cValue,1,1,m_pFile);
	specCurrent.m_fZ = (float)cValue;
	nread = fread((void *)&cValue,1,1,m_pFile);
	nread = fread((void *)pLine,1,(int)cValue,m_pFile);
	pLine[cValue] = '\0';
	specCurrent.m_strDescription = pLine;
	fValue = 0.0;
	nread = fread((void *)&fValue,sizeof(float),1,m_pFile);
	float fIntensity = fValue;
	cValue = 0;
	nread = fread((void *)&cValue,1,1,m_pFile);
	size_t tSize = (size_t)cValue;
	fValue = 0.0;
	nread = fread((void *)&fValue,sizeof(float),1,m_pFile);
	nread = fread((void *)&cValue,1,1,m_pFile);
	float fScale = fValue;
	size_t a = 0;
	nread = fread((void *)&sValue,2,1,m_pFile);
	iValue = (unsigned int)sValue;
	miCurrent.m_fM = (float)iValue/fScale;
	specCurrent.m_vMI.push_back(miCurrent);
	a++;
	while(a < tSize)	{
		nread = fread((void *)&sValue,2,1,m_pFile);
		iValue += (unsigned int)sValue;
		miCurrent.m_fM = (float)iValue/fScale;
		specCurrent.m_vMI.push_back(miCurrent);
		a++;
	}
	a = 0;
	double dSum = 0.0;
	char cMax = 0;
	while(a < tSize)	{
		nread = fread((void *)&cValue,1,1,m_pFile);
		specCurrent.m_vMI[a].m_fI = (float)cValue;
		if(cMax < cValue)	{
			cMax = cValue;
		}
		dSum += (double)cValue;
		a++;
	}
	delete pLine;
	specCurrent.m_vdStats.push_back(dSum*fIntensity);
	specCurrent.m_vdStats.push_back((double)cMax*fIntensity);
	specCurrent.m_vdStats.push_back((double)fIntensity);
	specCurrent.m_tId = m_tId;
/*
 * copy the temporary mspectrum into the input reference
 */
/*
 * if the file is finished, return false, otherwise return true 
 */
	if(feof(m_pFile))	{
		fclose(m_pFile);
		return false;
	}
	_m = specCurrent;
	return true;
}

bool loadcmn::open(string &_s)
{
	m_tId = 1;
/*
 * copy the file pathname into m_strPath
 */
	m_strPath = _s;
/*
 * try to open the file and bail out if it isn't available
 */
	m_pFile = fopen(m_strPath.c_str(),"rb");
	if(m_pFile == NULL)	{
		return false;
	}
/*
 * check the open file to test for the value of the end-of-line character
 */
	char *pLine = new char[m_tSize];
	size_t nread = fread((void *)pLine,1,256,m_pFile);
	pLine[255] = '\0';
	if(strstr(pLine,"CMN ") != pLine)	{
		fclose(m_pFile);
		return false;
	}
	delete pLine;
	return true;
}
// loaddta::open_force
// forces the opening of a DTA format file
bool loadcmn::open_force(string &_s)
{
	return open(_s);
}

#ifdef HAVE_PWIZ_MZML_LIB
//
// use the ProteoWizard library to load mzML and anything else it can handle
//
#undef __int64_t // conflict with boost
#undef uint32_t  // conflict with boost
#undef uint64_t  // conflict with boost
#include <exception>
#include <RAMPAdapter.hpp>
loadpwiz::loadpwiz(void):handler(NULL)
{
};
loadpwiz::~loadpwiz(void )
{
	delete handler;
};
 
bool loadpwiz::get(mspectrum &_m) {
	if (m_chargestack.size()) {
		// we have another guess at charge state
		_m = m_chargestack[m_chargestack.size()-1];
		m_chargestack.pop_back();
		return true;
	}
	bool result = false;
	// read the next MS2 spectrum
	while ((!result) && (this->m_tId < handler->scanCount())) {
		ScanHeaderStruct hdr;
		handler->getScanHeader(this->m_tId, hdr, false); // false = don't read peaks yet
		if (hdr.msLevel==2) {
			std::vector<double> peaks;
			handler->getScanPeaks(this->m_tId, peaks);
			/* load the parent ion charge and intensity	*/
			_m.m_fZ = (float)hdr.precursorCharge;
			_m.m_fI = (float)hdr.precursorIntensity;

			/* load the m/z and intensity values */
			_m.m_vMI.resize(peaks.size()/2);
			int nn=0;
			for (int n=0;n<(int)_m.m_vMI.size();n++) {
				_m.m_vMI[n].m_fM = (float)peaks[nn++];
				_m.m_vMI[n].m_fI = (float)peaks[nn++];
			}
			/* note retention time */
			if (hdr.retentionTime>0) {
				char buf[256];
				snprintf(buf,sizeof(buf),"%f",hdr.retentionTime);
				_m.m_strRt=buf;
			}

			const double dProton = 1.007276; // why isn't this a global constant?
			if (!_m.m_fZ) {  
				// undeclared charge state - if it's not obviously singly charged 
				// then submit as 2+ and 3+ a la mzXML parser, from which this
				// charge guessing code was lifted
				// guess the charge state of the precursor
				// from the ratio of the integrals of the intensities below and above m_precursorMz
				float intBelow = 0;
				float intTotal = 0;

				size_t length = _m.m_vMI.size();
				for(size_t n = 0 ; n < length ; n++)
				{
					intTotal += _m.m_vMI[n].m_fI;
					if(_m.m_vMI[n].m_fM < hdr.precursorMZ)
						intBelow += _m.m_vMI[n].m_fI;   
				}

				// There is no particular reason for the 0.95. It's there just
				// to compensate for the noise.... 
				if(intTotal == 0.0 || intBelow/intTotal > 0.95) {
					_m.m_fZ = 1;
				} else {
					// doesn't appear to be singly charged - emit as 2+ and 3+
					_m.m_fZ = 3;
					_m.m_dMH = (hdr.precursorMZ - dProton)*_m.m_fZ + dProton;
					_m.m_tId = hdr.acquisitionNum+100000000; // attempt at unique ID
					// next call here will get the 3+, this call will return the 2+
					m_chargestack.push_back(_m);
					_m.m_fZ = 2;
				}
			}
			_m.m_dMH = (hdr.precursorMZ - dProton)*_m.m_fZ + dProton;
			_m.m_tId = hdr.acquisitionNum; // use machine's scan number as ID

			result = true;
		}
		this->m_tId++;
	}
	return result;
};

#include <stdexcept>

bool loadpwiz::open(string &_s)
{
	try {
		handler = new pwiz::msdata::RAMPAdapter(_s);
	}
	catch (std::runtime_error &e) {
		// unsupported format, presumably
		std::string what(e.what());
		delete handler;
		handler = NULL;
	}
	if (handler) { // pwiz can read it
		this->m_strPath = _s;
		this->m_tId = 0;
	}
	return (NULL!=handler);
};
bool loadpwiz::open_force(string &_s)
{
	return open(_s);
};

#endif // ifdef HAVE_PWIZ_MZML_LIB
