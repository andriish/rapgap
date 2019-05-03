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
    //int getorig_(int &a);
    //int ncount_hepmc3 = 0;
    //char * outfile_hepmc3 = getenv ("HEPMC3OUT");

    //HepMC3::WriterAscii     ascii_io_hepmc3(outfile_hepmc3);
    //HepMC3::GenCrossSection cross_hepmc3;



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

HepMC::GenEvent*  gev;


  int rivetinit_() {
    if ( rivet ) return 0;
    rivet = new Rivet::AnalysisHandler(rname);
    rivet->setIgnoreBeams(igBeam);
    Rivet::addAnalysisLibPath(".");
    for(int i = 0, N = preloads.size(); i < N; ++i)
      rivet->readData(preloads[i]);
    for (set<string>::iterator it = analyses.begin();
      it != analyses.end(); ++it) {
      rivet->addAnalysis(*it);
    }
    rivet->init(*gev);
    return 0;
  }


  int rivetdone_() {
    if ( !rivet ) return 0;
    rivet->finalize();
    rivet->writeData(filename);
    delete rivet;
    rivet = 0;
    return 0;
  }



    void convrivet_(int & ievent, int & iproc, double & xsec, double & xsece,
                    int& flav1, int& flav2,double &  x1,double &  x2,double &  q2pdfeval,
                    double &  xf1mom, double & xf2mom ,int&  pdf1,int& pdf2
                   ) {

       
    }
}
