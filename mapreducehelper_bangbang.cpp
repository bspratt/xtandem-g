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

#ifdef HAVE_MULTINODE_TANDEM // support for Hadoop and/or MPI?

#include <iostream>
#undef READ_GZIP // messes with ifstream etc
#include "stdafx.h"
#include <algorithm>
#include <set>
#include "msequence.h"
#include "msequencecollection.h"
#include "msequenceserver.h"
#include "msequtilities.h"
#include "mspectrum.h"
#include "loadmspectrum.h"
#include "xmlparameter.h"
#include "mscore.h"
#include "mprocess.h"
#include "saxbiomlhandler.h"
#include "mbiomlreport.h"
#include "mrefine.h"
#include "serialize.h"
#include "timer.h"
#include "mapreducehelper.h"

#include <unistd.h>
#include <string>
#include <sstream>


//
// the contents of this file are largely lifted from 
// ownercompute.cxx in the X!!Tandem project.
//


#ifdef XBANGBANG
#include "mpi.h"
#else
void no_mpi() {
	cerr << "not compiled for X!!Tandem MPI use, quitting\n";
	exit(1);
}
#endif

// set up MPI usage, if we're built for that
void mapreducehelper::init_mode_bangbang(int argc,char **argv) {
#ifndef XBANGBANG
	no_mpi();
#else
	mode_bangbang_ = true;
	MPI::Init(argc,argv);
	num_processes_ = MPI::COMM_WORLD.Get_size();	
	ID_ = MPI::COMM_WORLD.Get_rank();
	if (!ID_) {
		t_ = new timer("master");
	}
#endif
}

/* This is necessary because of a conflict between C++ headers and mpi headers
   See the comment in mpicxx.h
*/
#ifdef XBANGBANG
#define GATHER 1
#define SCATTER 2
#endif

void mapreducehelper::abort() {
#ifdef XBANGBANG
	MPI::COMM_WORLD.Abort(1);
#endif
	exit(mode_mapreduce()?0:1);  // yes, an error but keep hadoop flow going so we can get logs easily
}

int mapreducehelper::mpi_rank() const { // rank of this process if in bangbang mode
#ifndef XBANGBANG
	no_mpi();
	return 0;
#else
    return ID_;
#endif
}
void mapreducehelper::mpi_gather(mprocess **pProcess) {
#ifndef XBANGBANG
	no_mpi();
#else

    if (mpi_rank() > 0) { /* workers */
      /* serialize my data */
      int len; // will be set to length of serialized data packet
      char *buf = serialize_process(*pProcess[mpi_rank()],&len,true);
      /* send to process 0 */
      MPI::COMM_WORLD.Send(buf, len, MPI::CHAR, 0, GATHER);
	  delete[] buf;
    } else { /* master */
      MPI::Status status;
      t_->addsplit("gather");
      for (int i=1; i<num_processes_; ++i) {
        t_->addsplit("probe");
        MPI::COMM_WORLD.Probe(MPI::ANY_SOURCE, GATHER, status);
        int source = status.Get_source();
        int len = status.Get_elements(MPI::CHAR);
        char *buf = new char[len];
        t_->addsplit("get data");
        MPI::COMM_WORLD.Recv(buf, len, MPI::CHAR, source, GATHER);
        t_->addsplit("got data");
        deserialize_process(buf, true, *pProcess[source]);
        delete[] buf;
      }
    }
#endif
}
void mapreducehelper::mpi_scatter(mprocess **pProcess) {
#ifndef XBANGBANG
	no_mpi();
#else
    char *buf;
    int len, i;
    MPI::Status status;
    std::string s;
    if (mpi_rank() > 0) { /* workers */
        MPI::COMM_WORLD.Probe(0, SCATTER, status);
        len = status.Get_elements(MPI::CHAR);
        buf = new char[len+1];
        MPI::COMM_WORLD.Recv(buf, len, MPI::CHAR, 0, SCATTER);
        deserialize_process(buf, true, *pProcess[mpi_rank()]);
        delete[] buf;
    } else { /* master */
      t_->addsplit("scatter");
      for (i=1; i<num_processes_; ++i) {
        /* serialize worker data */
	    int len;
cout << "serialize a process with address" << (long)pProcess[i] << " \n";
cout << "serialize a process with " << pProcess[i]->m_vSpectra.size() << " spectra\n";
        buf = serialize_process(*(pProcess[i]),&len,true);
        /* send to process i */
t_->addsplit("send");
        MPI::COMM_WORLD.Send(buf, len, MPI::CHAR, i, SCATTER);
t_->addsplit("ok,ish");
		delete[] buf;
      }
      t_->addsplit("send data");
    }
#endif
}
void mapreducehelper::mpi_finalize() {
#ifndef XBANGBANG
	no_mpi();
#else
	int rank = mpi_rank();
	MPI::Finalize();
	if (rank > 0) {
		exit(0);
	}
#endif
}

#endif // HAVE_MULTINODE_TANDEM

