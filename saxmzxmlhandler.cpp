/************************************************************
 * SAXMzxmlHandler.cpp
 *
 * Premiere version janvier 2005
 * Patrick Lacasse
 * placasse@mat.ulaval.ca
 *
 * 3/11/2005 (Brendan MacLean): Use eXpat SAX parser, and create SAXSpectraHandler
 *
 * See http://sashimi.sourceforge.net/software_glossolalia.html for
 * mzXML schema information.
 *
 * Inspired by DtaSAX2Handler.cpp
 * copyright            : (C) 2002 by Pedrioli Patrick, ISB, Proteomics
 * email                : ppatrick@systemsbiology.org
 * Artistic License granted 3/11/2005
 *******************************************************/

#include "stdafx.h"
#include "saxmzxmlhandler.h"

SAXMzxmlHandler::SAXMzxmlHandler( vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m)
	: SAXSpectraHandler(_vS, _sC, _m)
{
	m_bInMsLevel2 = false;
	m_bInPrecursorMz =  false;
	m_bInPeaks = false;
	//added this (true by default)
	m_bNetworkData = false;
}

SAXMzxmlHandler::~SAXMzxmlHandler()
{
}

void SAXMzxmlHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
	if(isElement("scan", el))
	{
		if((m_cidLevel = atoi(getAttrValue("msLevel", attr))) == 2)
		{
			m_bInMsLevel2 = true;

			reset();	// Clean up for the next scan

			m_scanNum = atoi(getAttrValue("num", attr));
			m_tId = m_scanNum;
			while(m_sId.find(m_tId) != m_sId.end())	{
				m_tId++;
			}
			m_sId.insert(m_tId);
			m_peaksCount = atoi(getAttrValue("peaksCount", attr));
			m_strRt = getAttrValue("retentionTime", attr);
		}
	}
	else if(isElement("peaks", el))
	{
		m_bInPeaks = true;
		m_bLowPrecision = (strcmp("64", getAttrValue("precision", attr)) != 0);
        m_bCompressed = false;
      const char *compressionType = getAttrValue("compressionType", attr);
      if ((*compressionType != '\0')&&strcmp(compressionType, "none"))
      {
          // TODO: Cause parse error with line number.
          if (strcmp(compressionType, "zlib") != 0)
          {
              cerr << "Unsupported compression type '" << compressionType << "'.\n";
              exit(EXIT_FAILURE);
          }
          m_bCompressed = true;

          const char *compressedLen = getAttrValue("compressedLen", attr);
          if (*compressedLen == '\0')
          {
              cerr << "Missing compressedLen attribute.\n";
              exit(EXIT_FAILURE);
          }

          m_lenCompressed = atoi(compressedLen);
      }
	}
	else if(isElement("precursorMz", el))
	{
		if (m_cidLevel < 3) {
			//don't read precursor data if ms level >= 3
		m_bInPrecursorMz = true;
		m_precursorCharge = atoi(getAttrValue("precursorCharge", attr));

		// test for the mzXML 3.1 possibleCharges attribute
		// ex: "7,11,13"
		string possibleCharges = getAttrValue("possibleCharges", attr);
		if (possibleCharges != "") {
		  // parse the comma-separated list of additional possible precursor charges, as may come from ETD data.
		  string::size_type token_begin = 0;
		  string::size_type token_end = string::npos;
		  bool done=false;
		  while (!done) {
		    token_end=possibleCharges.find_first_of(',', token_begin);
		    string charge;
		    if (token_end == string::npos) {
		      charge=possibleCharges.substr(token_begin);
		      done = true;
		    }
		    else {
		      charge=possibleCharges.substr(token_begin, token_end-token_begin);
		      token_begin = token_end+1;
		    }
		    if (charge.size()>0) {
		      m_viPossiblePrecursorCharges.push_back(atoi(charge.c_str()));
		    }
		  }
		}



	}
	}
}

void SAXMzxmlHandler::endElement(const XML_Char *el)
{
	if(isElement("peaks", el))
	{
		processData();
		m_bInPeaks = false;
	}
	else if(isElement("precursorMz", el))
	{
		processData();
		m_bInPrecursorMz = false;
	}
	else if(isElement("scan", el) && m_bInMsLevel2 == true)
	{
		// only add a spectrum without charge (which will lead
		// to internal xtandem charge state guessing) if there
		// were no values parsed from *both* "precursorCharge"
		// or "possibleCharges"
		if ( (m_precursorCharge == 0) && (m_viPossiblePrecursorCharges.size() == 0) ) {
		  // add spectrum, with precursorMz charge
		  pushSpectrum();
		}
		
		else {

		  // add the spectrum with the m_precursorCharge value
		  int originalPrecursorMZ = m_precursorCharge; // do other pushSpectrum calls change this?
		  pushSpectrum(m_precursorCharge);
		  
		  // are there any multiple precursor charges from mzXML 3.1's
		  // possibleCharges?
		  if (m_viPossiblePrecursorCharges.size() > 0) {
		    size_t originalId = m_tId;
		    for (vector<int>::iterator i = m_viPossiblePrecursorCharges.begin();
			 i != m_viPossiblePrecursorCharges.end();
			 ++i) {
		      int z = *i;
		      if (z != originalPrecursorMZ) { // no need to duplicate if already added
			m_tId += 100000000;
			pushSpectrum(z);
		      }
		    }
		    m_tId = originalId;
		  }
		}

		m_bInMsLevel2 = false;
	}
}

void SAXMzxmlHandler::characters(const XML_Char *s, int len)
{
	if ((m_bInPeaks && m_cidLevel == 2) ||
		(m_bInPrecursorMz))
	{
		m_strData.append(s, len);
	}
}

void SAXMzxmlHandler::processData()
{
	if( m_bInPeaks && m_cidLevel == 2)
	{
		pushPeaks();
	}
	else if (m_bInPrecursorMz)
	{
		if (m_cidLevel < 3) {
		m_precursorMz = atof(m_strData.data());
		}
	}

	m_strData.clear();
}
