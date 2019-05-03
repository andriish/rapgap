#include <iostream>
#include "HepMC3/WriterAscii.h"
#include "HepMC3/GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include <set>
#include <vector>
#include <string>
#include "Rivet/Rivet.hh"
using namespace std;

extern "C" {
extern HepMC::GenEvent* evt;
  /**
   * The name of the file where the histograms are dumped.
   */
   string filename;

  /**
   * Analyses with optional analysis parameters.
   */
  set<string> analyses;

  /**
   * The names of YODA files to preload.
   */
  vector<string> preloads;

  /**
   * The Rivet object.
   */
  Rivet::AnalysisHandler * rivet;


  /**
   * The Rivet run name
   */
  string rname;

  /**
   * Ignore beams flag.
   */
  bool igBeam=true;


  int rivetinit_(char* rname1) {
    if ( rivet ) return 0;
    rname=string(rname1);
    rivet = new Rivet::AnalysisHandler(rname);
    rivet->setIgnoreBeams(igBeam);
    Rivet::addAnalysisLibPath(".");
    for(int i = 0, N = preloads.size(); i < N; ++i)
      rivet->readData(preloads[i]);
    for (set<string>::iterator it = analyses.begin();
      it != analyses.end(); ++it) {
      rivet->addAnalysis(*it);
    }
    
    return 0;
  }
int rivetinitfirstevent_()
{
rivet->init(*evt);
}
int rivetrun_(){
    
    rivet->analyze(*evt);
  }


  int rivetadd_(char* ana)
  {
	 analyses.insert(std::string(ana));
	 return   analyses.size();
  }  
  int rivetdone_(char* filename1) {
    if ( !rivet ) return 0;
    filename=std::string(filename1);
    rivet->finalize();
    printf("->%s<-\n",filename.c_str());
    rivet->writeData(filename);
    delete rivet;
    rivet = 0;
    return 0;
  }
}
