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

#ifndef MAPREDUCEHELPER_INCL___
#define MAPREDUCEHELPER_INCL___

class timer;

enum eMapReduceMode {NONE=0, MAPPER_1, REDUCER_1, MAPPER_2, REDUCER_2, MAPPER_3, REDUCER_3};

void mapreduce_heartbeat_until_exit(); // just emit heartbeat on stderr till thread calls exit

class mapreducehelper {
private:
  int num_processes_; // populated by reading stdin, lines of format <int process>,<int total_processes>
  vector<int> process_ids_; // populated by reading stdin, lines of format <int process>,<int total_processes>
  eMapReduceMode mode_mapreduce_; // NONE, MAPPER, REDUCER 
  int step_mapreduce_; // which map/reduce step
  std::string modetext_; // "mapper1", "reducer3" etc
  int ID_; // unique ID for this mapper or reducer (read from input stream), or MPI rank
  std::string ID_str_; // string representation of ID for mapper or reducer
  int n_processes_written_; // how many processes have been serialized
  int requested_mapper_count_; // how many mappers we actually hoped for
  std::vector<int> task_mapreduce_; // for local debugging, can specify MAPPER_1.2.3 to mean "pretend you are task 2 and 3 for MAPPER_1 step"
  timer *t_;
  
  // for mapreduce implementation we redirect std::cout to std::cerr,
  // and hook stdout up to our own std::ostream mapreduce_out
  ostream *mapreduce_out_;

  // for large data sets Hadoop can fail trying to sort data coming in on stdio.
  // so instead we write to HDFS and refer to those files in stdio
  std::string outdir_;

  // URL for final file (in AWS, normally s3n://something)
  std::string report_URL_;

  // for X!!Tandem implementation
  bool mode_bangbang_;

public:
  mapreducehelper(int &argc, char **&argv);
  ~mapreducehelper();
  // get processes, if any, from stdin (or a stdin reference to a file)
  int read_processes(std::vector<mprocess *> &result); // return number of processes read
  
  void send_processes_to_next_task(mprocess **pProcessList,int processCount); 

  void copy_report_to_url(mprocess *p); // get the filename from process

  eMapReduceMode mode_mapreduce() const; // return nonzero if mode_mapreduce_ is mapper or reducer

  static std::string timestamp(); // can be useful for debugging

  bool mode_bangbang() const; // return true iff running in MPI X!!Tandem mode
  int mode_bangbang_client() const; // return nonzero iff a client running in MPI X!!Tandem mode
#ifdef _DEBUG
  void exit_mode_mapreduce(); // for debugging
#endif
  bool accept_this_mapreduce_input_line(const std::string &line); // return TRUE if the indicated task number is ours (a local debug thing)
  void report();
  void abort(); // exit with nonzero code, kill MPI if in bangbang mode

  static void sort_process_spectra(mprocess &process); // gives consistent order to mimic single thread operation

  std::string numprocs() const; // formatted number of known mapreduce processes, or size of MPI ring
  int num_processes() const {
	  return num_processes_;
  }

  // X!!Tandem specific stuff
  // set up MPI usage, if we're built for that
  void init_mode_bangbang(int argc,char **argv);
  int mpi_rank() const; // rank of this process if in bangbang mode
  void mpi_gather(mprocess **pProcess);
  void mpi_scatter(mprocess **pProcess);
  void mpi_finalize();

  // hacked-up version of tandem main routine for mapreduce and X!!tandem
  int mapreducehandler(int argc, char* argv[]); // in mapreducehandler.cpp

};

mapreducehelper &the_mapreducehelper(); // this is a singleton


#endif // MAPREDUCEHELPER_INCL___

#endif // HAVE_MULTINODE_TANDEM 
