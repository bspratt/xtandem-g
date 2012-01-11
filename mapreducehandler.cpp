/*

 MapReduce implementation of X!Tandem
 Copyright (C) 2010 Insilicos LLC All Rights Reserved

 Also includes MPI implementation based on X!!Tandem which is
 Copyright (C) 2007 Robert D Bjornson, all rights reserved

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

// File version: 2004-02-01
// File version: 2004-03-01

/*
   Modified 2010 Brian Pratt for X!!Tandem MapReduce version.
*/

/*
	tandem.cpp provides a command line interface to the main processing class, mprocess.
	A single command line parameter is accepted, which is the file name for an XML class
	containing the parameters for performing the protein modelling session. The input file 
	also contains the name of an output file, which will contain the output from the session.
*/

#ifdef HAVE_MULTINODE_TANDEM  // support for Hadoop and/or MPI?

#include "stdafx.h"
#include <sys/timeb.h>
#include <ctime>
#include <sstream>
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"
#include "timer.h"
#include "mapreducehelper.h"
#include "serialize.h"

static int MAXTHREADS=16;

#ifndef X_P3

/*
 * windows.h and the definition of ProcessThread are necessary for multithreading
 * in the Windows 32 environment. UNIX platforms that use POSIX threads use an
 * alternate version of this file.
 */
#ifdef MSVC
#include "windows.h"
#define THREAD_RETURN_TYPE_ DWORD 
#define THREAD_API_ WINAPI 
#define THREAD_VOIDPOINTER_ LPVOID
#else
#include <pthread.h>
#define THREAD_RETURN_TYPE_ void * 
#define THREAD_API_ 
#define THREAD_VOIDPOINTER_ void *
#endif
#define THREAD_DECL_ THREAD_RETURN_TYPE_ THREAD_API_ 
THREAD_DECL_ MRProcessThread(THREAD_VOIDPOINTER_ pParam);
THREAD_DECL_ MRRefineThread(THREAD_VOIDPOINTER_ pParam);
THREAD_DECL_ oldSkoolHeartbeatThread(THREAD_VOIDPOINTER_ pParam);

//
// routine to tidy up paths we flattened for use with Hadoop cachefile systems
//
static void tidy_path_string(std::string &str) {
	size_t pos;
	while ((pos=str.find("__sl__"))!=std::string::npos) {
		str.replace(pos,6,1,'/');
	}
	while ((pos=str.find("__cln__"))!=std::string::npos) {
		str.replace(pos,7,1,':');
	}
	// look for cygwin "/cygdrive/x/" weirdness, go back to straight windows "x:/"
	while ((pos=str.find("/cygdrive/"))!=std::string::npos) {
		std::string cygdrive = str.substr(pos,12);
		std::string drive = cygdrive.substr(10,1)+":/";
		str.replace(pos,12,drive);
	}
}

int mapreducehelper::mapreducehandler(int argc, char* argv[])
{
	// mapreduce steps are these (they correspond to traditional thread model):
	// mapper1: 
	//	input=one or more lines of text
	//  action=eat input text just to show you're alive
	//	output=single line of text (key = 1) just to show you're alive 
	// reducer1: (single reducer)
	//	input=one line per active mapper
	//  action=load, condition and load-balance spectra across counted mappers
	//  output=serialized spectra-count-rebalanced process objects for multiple mapper2's to process
	// mapper2:
	//	input=one or more serialized process objects
	//  action=run the "process" thread
	//	output=serialized process object with processed spectra, all mappers emit same key (single reducer)
	// reducer2: (single reducer)
	//	input=process objects with processed spectra
	//  action=make the processes aware of each other's shared state if refining, else stitch together final result
	//  output=serialized process objects for mapper3 to process if refining, else final result
	// mapper3:
	//	input=process objects with processed spectra
	//  action=refinement
	//	output=serialized process object with processed spectra, all with same key
	// reducer3: (single reducer)
	//  action: stitch together the results
	//	output=final result

	bool verbose = true;
	if (mode_mapreduce()) { 
		MAXTHREADS=1; // multiprocessor, not multithread
	} else if (mode_bangbang_client()) {  
		verbose = false;
		MAXTHREADS = num_processes()+1; // X!!Tandem client node
	} else {  // X!!Tandem head node
		if (mode_bangbang()) {  
			MAXTHREADS = num_processes()+1; // X!!Tandem head node
		} else {
			// this isn't happening
			cerr << "internal error!" ;
			exit(-2);
		}
		if(argc < 2 || argc > 1 && strstr(argv[1],"-L") == argv[1] || argc > 1 && strstr(argv[1],"-h") == argv[1])	{
			return -1;
		}
	} // end if X!!Tandem head node, or basic operation
	cout << timestamp() << "X! TANDEM " << VERSION << "\n\n";
	unsigned long lMaxThreads = MAXTHREADS;

	if (REDUCER_1 == mode_mapreduce()) {
		// count the number of output lines from MAPPER_1
		lMaxThreads = 0;
		std::cerr <<  timestamp() << "read lines from mapper_1 to estabish mapper count\n";
		std::string line;
		std::set<std::string> unique; // each mapper reports its hostname
		while (getline(cin,line)) {
			lMaxThreads++;
			unique.insert(line); // identify the number of unique hostnames
		}
		if (!lMaxThreads) {
			std::cerr << "no mappers?!?\n" ;
			exit(0); // yes, an error, but keep the flow going so we can retrieve logs
		}
		std::cerr <<  timestamp() << lMaxThreads << " mappers reporting in with " << unique.size() << " unique hostnames\n";
		// hadoop tends to give more mappers than we asked for in phase 1
		// be conservative and only plan on the ones we requested in phase 2
		// what we're really trying to do is avoid expecting more than we are given
		if ((requested_mapper_count_ > 0) && (requested_mapper_count_ < (int)lMaxThreads)) {
			lMaxThreads = requested_mapper_count_;
		}
		// create MR_DUP_FACTOR times as many
		// output files as nodes so that failure 
		// of any one node doesn't double the load of another
		if (lMaxThreads > 1) { // no point in single node case
			const int MR_DUP_FACTOR=4; //6;
			lMaxThreads *= MR_DUP_FACTOR;
		}
	}

	/*
	* Create an mprocess object array
	*/
	mprocess **pProcess = new mprocess*[lMaxThreads];
	if(pProcess == NULL)	{
		cout << timestamp() << "An error was detected creating the processing objects.\nPlease contact a GPM administrator.\n";
		if (mode_bangbang()) {  // X!!Tandem MPI mode
			return -2;
		}
		exit(0); // even though this an error, for hadoop we want to exit gracefully so we can pass error text back to user
	}
	unsigned long a = 0;
	while(a < lMaxThreads)	{
		pProcess[a] = NULL;
		a++;
	}
    /* moved these declarations here so that they are seen whether in mapreduce mode or not */
    int x=0;
	void *vp=NULL;
#ifdef MSVC
	DWORD *pId = new DWORD[lMaxThreads];
	DWORD dCount;
	DWORD dwTime=0;
	int iTics = 0;
#else
	int *pId = new int[lMaxThreads];
	pthread_t pThreads[lMaxThreads];
	int	dCount;
#endif
	pProcess[0] = new mprocess;
	if (!mode_mapreduce() || (REDUCER_1==mode_mapreduce())) { // need to read raw spectra?
		if (verbose) { 
			cout << timestamp() << "Loading spectra .";
			cout.flush();
		}
		/*
		* Initialize the first mprocess object with the input file name.
		*/
		if(!pProcess[0]->load(argv[1]))	{
			cout << "\n\n" << timestamp() << "An error was detected while loading the input parameters.\nPlease follow the advice above or contact a GPM administrator to help you.";
			delete pProcess[0];
			delete pProcess;
			if (mode_bangbang()) {  // X!!Tandem MPI mode
				return -4;
			}
			exit(0); // even though this an error, for hadoop we want to exit gracefully so we can pass error text back to user
		}
		pProcess[0]->set_thread(0);
		if (mode_bangbang()) {  // X!!Tandem MPI mode
			string strKey="spectrum, threads";
			string strValue=numprocs();
			cout << timestamp() << numprocs() << " workers ready\n";
			pProcess[0]->m_xmlValues.set(strKey, strValue);
			pProcess[0]->set_threads(num_processes());
		}  else {
			pProcess[0]->set_threads(lMaxThreads);
		}
		if (verbose) {  
			cout << timestamp() << "loaded.\n";
		}
	} // end if need to read raw spectra

	if (mode_bangbang() || (REDUCER_1 == mode_mapreduce())) {

		if(pProcess[0]->m_vSpectra.size() == 0)	{
			if (verbose) {  
				cout << timestamp() << "No input spectra met the acceptance criteria.\n";
				cout.flush();
			}
			if (mode_bangbang()) { // for mapreduce we'd like a clean exit even if empty - there may be other jobs to follow
				delete pProcess[0];
				delete pProcess;
				return 1;
			}
		}
		if (verbose) {  
			cout << timestamp() << "Spectra matching criteria = " << (unsigned long)pProcess[0]->m_vSpectra.size() << "\n";
			cout.flush();
#ifdef PLUGGABLE_SCORING
			cout << timestamp() << "Pluggable scoring enabled.\n";
#endif
		}
		/*

		* Start the mprocess object and wait for it to return.

		*/
		unsigned long lThread =	pProcess[0]->get_thread();
		unsigned long lThreads = pProcess[0]->get_threads();

		if (mode_mapreduce()) {  // check for mapper count too high for spectra count
			const int min_spectra_per_thread = 100; // this value is completely arbitrary
			if ((pProcess[0]->m_vSpectra.size() / lThreads) < min_spectra_per_thread) {
				lMaxThreads = (pProcess[0]->m_vSpectra.size() / min_spectra_per_thread);
				if (!lMaxThreads) {
					lMaxThreads = 1;
				}
			}
		}

		if(lThreads	> lMaxThreads)	{
			if (mode_bangbang_client()==1) { /* Parallel version must have same # threads as processes, so this X!Tandem fix won't work for us. */
				cerr << "ERROR: MPI allocation of " << lThreads << " exceeded maximum number of threads (" << lMaxThreads << ").\n";
				cerr << "Change MAXTHREADS in tandem.cpp, and recompile.\n";
				abort();
			}
			lThreads = lMaxThreads;
			pProcess[0]->set_threads(lThreads);
		}
		if(pProcess[0]->m_vSpectra.size() <	lThreads)	{
			if (mode_bangbang_client()==1) { /* Parallel version can't handle this case either */
				cerr << "ERROR: MPI allocation of " << lThreads << " exceeded number of spectra (" << pProcess[0]->m_vSpectra.size() << ").\n";
				cerr << "Rerun with fewer mpi processes.\n";
				abort();
			}
			lThreads = (unsigned long)pProcess[0]->m_vSpectra.size();
			if(lThreads	< 1)		{
				lThreads = 1;
			}
			pProcess[0]->set_threads(lThreads);
		}
		dCount = lThreads -	1;
		long lSpectra =	lThreads + (long)pProcess[0]->m_vSpectra.size()/lThreads;
		bool bSpectra =	true;
		if (verbose) {  
			cout << timestamp() <<	"Starting threads .";
			cout.flush();
		}
		if(lThread != 0xFFFFFFFF)		{
			while(dCount > 0)	{
				pProcess[dCount] = new mprocess;
				pProcess[dCount]->set_thread(dCount);						 
				/*

				* initialize the new mprocess objects with	the	spectra	already	loaded into	the	first mprocess

				*/
				pProcess[dCount]->m_vSpectra.reserve(lSpectra);
				dCount--;
			}
			size_t tCount = pProcess[0]->m_vSpectra.size();
			size_t tProcesses = lThreads;
			size_t tRing = 0;
			vector<mspectrum> vZero;
			vZero.reserve(lSpectra);
			do	{
				if(tRing == 0)	{
					vZero.push_back(pProcess[0]->m_vSpectra.back());
				}
				else	{
					pProcess[tRing]->m_vSpectra.push_back(pProcess[0]->m_vSpectra.back());
				}
				tRing++;
				pProcess[0]->m_vSpectra.pop_back();
				if(tRing == tProcesses)	{
					tRing = 0;

				}
			}	while(pProcess[0]->m_vSpectra.size() != 0);
			pProcess[0]->m_vSpectra.reserve(vZero.size());
			do	{
				pProcess[0]->m_vSpectra.push_back(vZero.back());
				vZero.pop_back();
			}	while(vZero.size() != 0);

			dCount = lThreads - 1;
			while(dCount > 0)		{
				if(!pProcess[dCount]->load(argv[1],pProcess[0]))	{
					cout <<	timestamp() <<"error pProcess->LoadParameters	returned error (main)\r\n";
					delete pProcess;
					if (mode_bangbang()) {  // X!!Tandem MPI mode
						return -4;
					}
					exit(0); // for hadoop we want to exit gracefully so we can pass error text back to user
				}
				if (mode_bangbang()) {  // X!!Tandem MPI mode
					string strKey="spectrum, threads";
					string strValue=numprocs();
					pProcess[dCount]->m_xmlValues.set(strKey, strValue);
					pProcess[dCount]->set_threads(num_processes());
				}
				dCount--;
				if (verbose) {  
					cout <<	".";
					cout.flush();
				}
			}
		}
		if (mode_mapreduce()) { // REDUCER_1
			cout << timestamp() << "send " << lThreads << " processes to mapper_2\n";
			// serialize the processes out to mapper_2 tasks
			send_processes_to_next_task(pProcess,lThreads);
			delete[] pProcess;
			cout << timestamp() << "done\n";
			exit(0);
		} // end if REDUCER_1

	} // end if not mapreduce, or in REDUCER_1

	switch (mode_mapreduce()) { // do we need to read in serialized processes?
		case MAPPER_2:
		case REDUCER_2: 
		case MAPPER_3:
		case REDUCER_3: 
			{
			// map reduce mode - read in serialized processes
			std::vector<mprocess *>processes;
			lMaxThreads = 0;
			std::cerr << timestamp() << "loading processes...";
			lMaxThreads+=read_processes(processes);
			std::cerr << "done.\n";
			delete[] pProcess;
			if (!lMaxThreads) {
				std::cerr << "no data\n";
				exit(0); // not necessarily an error in mapreduce context
			}
			pProcess = new mprocess*[lMaxThreads+1];
			std::set<std::string> unique; // each mapper reports its hostname
			for (int n=lMaxThreads;n--;) {
				pProcess[n] = processes[n];
				unique.insert(pProcess[n]->m_hostname); // identify the number of unique hostnames
			}
			pProcess[lMaxThreads] = NULL; // null terminate the list
			if (lMaxThreads > 1) {
				std::cerr <<  timestamp() << lMaxThreads << " processes loaded from " << unique.size() << " unique hostnames\n";
			}
			dCount=lMaxThreads;
		} break;
		default:
			break;
	} // end do we need to read in serialized processes?

	// time to process spectra?
	if (MAPPER_2 == mode_mapreduce()) { 
		if (lMaxThreads > 1) {
			cout << timestamp() << "Combining processes:\n";
			cout.flush();
			size_t nspecs = 0;
			for (size_t m=lMaxThreads;m--;) {
				nspecs += pProcess[m]->m_vSpectra.size();
			}
			pProcess[0]->m_vSpectra.reserve(nspecs);
			for (size_t n=lMaxThreads;n-->1;) { // in case we get handed multiple jobs
				// combine spectra into single process
				while (pProcess[n]->m_vSpectra.size()) {
					pProcess[0]->m_vSpectra.push_back(pProcess[n]->m_vSpectra.back());
					pProcess[n]->m_vSpectra.pop_back();
				}
				delete pProcess[n];
				pProcess[n] = NULL;
			}
			pProcess[0]->m_tSpectraTotal = pProcess[0]->m_vSpectra.size();
		}
		cout << timestamp() << "Computing models on " << pProcess[0]->m_vSpectra.size() << " spectra:\n";
#ifdef MSVC
		CreateThread(NULL,0,MRProcessThread,(void *)pProcess,0,&pId[0]);
#else
		pthread_create(&pThreads[0],NULL,MRProcessThread,(void*)pProcess);
#endif
		mapreduce_heartbeat_until_exit(); // hadoop seems happier when heartbeat is in main thread and processing is in child
  } else if (mode_bangbang()) { // X!!Tandem

	  dCount = 1;
	  /*
	  * Initialize more mprocess objects, if lThread is not 0xFFFFFFFF, which signifies default single
	  * threaded operation.
	  */
	  unsigned long lThread =	pProcess[0]->get_thread();
	  unsigned long lThreads = pProcess[0]->get_threads();

	  if (verbose) {
		  cout << timestamp() << "push data to workers... ";
	  }
	  mpi_scatter(pProcess); // push master's data out to workers
	  if (verbose) {
		  cout << timestamp() << "process... ";
	  }
	  pProcess[mpi_rank()]->process();
	  if (verbose) {
		  cout << timestamp() << "gather results... ";
	  }
	  mpi_gather(pProcess);
	  /* RDB this happens as a side effect in the non-OC code, so I'll do it here */
	  dCount=lThreads;
  } // end if time to process spectra
  if (mode_bangbang() || (REDUCER_2==mode_mapreduce())) { // do we need to mingle all processes?

	  pProcess[0]->merge_spectra();
	  a = 1;
	  /*
	  * merge the results into the first object
	  */
	  while(a < dCount)	{
		  pProcess[0]->merge_map(pProcess[a]->m_mapSequences);
		  pProcess[0]->merge_spectra(pProcess[a]->m_vSpectra);
		  a++;
	  }
	  a = 1;
	  pProcess[0]->load_sequences();
	  while(a < dCount)	{
		  pProcess[a]->merge_map(pProcess[0]->m_mapSequences);
		  pProcess[a]->m_vseqBest = pProcess[0]->m_vseqBest;
		  a++;
	  }
  } // end do we need to mingle all processes

	bool bRefine = true;
	if (mode_mapreduce() && (MAPPER_3 != mode_mapreduce()) && (REDUCER_3 != mode_mapreduce())) { 
		// mapper 3 does refinement - do we need that?
		std::string strKey("refine");
		std::string strValue;
		pProcess[0]->m_xmlValues.get(strKey,strValue);
		if (strValue == "yes") {
			send_processes_to_next_task(pProcess,dCount);
			exit(0); // done with REDUCER_2
		} else {
			bRefine = false;  // we can finish up here
		}
	}

	// time for refinement?
	if (!bRefine) {
		; // do nothing
	} else if (MAPPER_3 == mode_mapreduce()) { 
		cout << timestamp() << "Model refinement on " << pProcess[0]->m_vSpectra.size() << " spectra:\n";
		cout.flush();
#ifdef MSVC
		CreateThread(NULL,0,MRRefineThread,(void *)pProcess,0,&pId[0]);
#else
		pthread_create(&pThreads[0],NULL,MRRefineThread,(void*)pProcess);
#endif
		mapreduce_heartbeat_until_exit(); // hadoop seems happier when heartbeat is in main thread and processing is in child
	} else 	if (mode_bangbang()) { 

		/*
		* Report the contents of the mprocess objects into an XML file as described
		* in the input file.
		*/
		if (verbose) {  
			cout << timestamp() << "Model refinement:\n";
			cout.flush();
		}
		dCount = 0;
		unsigned long lThread =	pProcess[0]->get_thread();
		unsigned long lThreads = pProcess[0]->get_threads();

		mpi_scatter(pProcess); // push master's data out to workers
		pProcess[mpi_rank()]->refine();
		mpi_gather(pProcess); // get workers' refinement data back to master
		mpi_finalize(); // workers exit here
		/* this happens as a side effect in the non-OC code, so I'll do it here */
		dCount=lThreads;
	} // end if time for refinement

	// if we got this far, we're either in X!! mode, or final mapreduce step reducer
	a = 1;
	/*
	* merge the results into the first object
	*/
	if(dCount > 1)	{
		cout << timestamp() << "Merging results:\n";
		cout.flush();
	}
	while(a < dCount)	{
		if(a == 1)	{
			cout << "\tfrom " << a+1;
		}
		else	{
			cout << a+1 << " ";
		}
		cout.flush();
		if(!pProcess[0]->add_spectra(pProcess[a]->m_vSpectra))	{
			cout << timestamp() << "adding spectra failed.\n";
		}
		pProcess[0]->merge_statistics(pProcess[a]);
		pProcess[a]->clear();
		pProcess[a]->m_mapSequences.clear();
		a++;
	}
	if(dCount > 1)	{
		cout << "\n\n";
		cout.flush();
	}
	cout.flush();
	cout << timestamp() << "Creating report:\n";
	cout.flush();
	// tidy up any path references for use back on calling system
	for (vector<mspectrum>::iterator iter = pProcess[0]->m_vSpectra.begin();
		         iter != pProcess[0]->m_vSpectra.end(); iter++) {
	    tidy_path_string(iter->m_strDescription);
		for (vector<msequence>::iterator s = iter->m_vseqBest.begin();
		    	s!=iter->m_vseqBest.end();s++) {
			tidy_path_string(s->m_strPath);
		}
	}
	for (vector<msequence>::iterator s = pProcess[0]->m_vseqBest.begin();
	    	s!=pProcess[0]->m_vseqBest.end();s++) {
		tidy_path_string(s->m_strPath);
	}
	// now write report
	pProcess[0]->report();
	size_t tValid = pProcess[0]->get_valid();
	size_t tUnique = pProcess[0]->get_unique();
	double dE = pProcess[0]->get_threshold();
	unsigned long lE = (unsigned long)(0.5+(double)tUnique/(1.0+1.0/dE));
	unsigned long lEe = (unsigned long)(0.5 + sqrt((double)lE));
	if(lEe == 0)	{
		lEe = 1;
	}
	if(dE <= 0.0)	{
		dE = 1.0;
	}
	cout << timestamp(); 
	cout << "\nValid models = " << (unsigned long)tValid << "\n";
	if(tUnique > 0)	{
		cout << "Unique models = " << (unsigned long)tUnique << "\n";
		cout << "Estimated false positives = " << lE << " &#177; ";
		cout << lEe << "\n";
	}
	lE = pProcess[0]->get_reversed();
	if(lE != -1)	{
		cout << "False positive rate (reversed sequences) = " << lE << "\n";
	}
	cout << "\n\n";

	// get the result off this local node and onto S3 where we can find it
	copy_report_to_url(pProcess[0]);  // no effect in non-mapreduce implementations

	/*
	* Delete the mprocess objects and exit
	*/
	a = 0;
	while(a < lMaxThreads)	{
		if(pProcess[a] != NULL)	{
			delete pProcess[a];
		}
		a++;
	}
	delete pProcess;
	delete pId;
	if (getenv("XTandem_TIMES"))
	    report();
	exit(0);
}
/*
 * Process thread is used to create the worker threads for each mprocess object
 */

THREAD_DECL_ MRProcessThread(THREAD_VOIDPOINTER_ _p)
{
	mprocess** pProcess = (mprocess**)_p; // a null terminated array of process pointers
	int n=0;
	while (*pProcess) { // in case we're handed more than one
		if ((*pProcess)->m_vSpectra.size()) {
			(*pProcess)->process();
		}
		pProcess++;
		n++;
	}
	the_mapreducehelper().send_processes_to_next_task((mprocess**)_p,n); // send to reducer, task #1
	exit(0);
	return (THREAD_RETURN_TYPE_)0;
}
THREAD_DECL_ MRRefineThread(THREAD_VOIDPOINTER_ _p)
{
	mprocess** pProcess = (mprocess**)_p; // a null terminated array of process pointers
	int n=0;
	while (*pProcess) { // in case we're handed more than one
		(*pProcess)->refine();
		pProcess++;
		n++;
	}
	the_mapreducehelper().send_processes_to_next_task((mprocess**)_p,n); // send to reducer, task #1
	exit(0);
	return (THREAD_RETURN_TYPE_)0;
}

#endif // end if not P3

#endif // HAVE_MULTINODE_TANDEM
