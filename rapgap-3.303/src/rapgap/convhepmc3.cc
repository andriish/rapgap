#include <iostream>
#include "HepMC3/WriterAscii.h"
#include "HepMC3/GenEvent.h"
#include "PythiaHelper.h"

using namespace std;

extern "C" {
    int getorig_(int &a);
    int ncount_hepmc3 = 0;
    char * outfile_hepmc3 = getenv ("HEPMC3OUT");

    HepMC3::WriterAscii     ascii_io_hepmc3(outfile_hepmc3);
    HepMC3::GenCrossSection cross_hepmc3;

    void convhepmc3_(int & ievent, int & iproc, double & xsec, double & xsece,
                    int& flav1, int& flav2,double &  x1,double &  x2,double &  q2pdfeval,
                    double &  xf1mom, double & xf2mom ,int&  pdf1,int& pdf2
                   ) {

        char * pPath;
        pPath = getenv ("HEPMC3OUT");
        if ( ncount_hepmc3 < 10) {
            if (pPath!=NULL) {
                cout << " env variable = " <<  pPath << endl;
            }
            else {
                cout << " NO HEPMC3OUT environment varibale set " <<endl;
                return ;
            }
            ++ncount_hepmc3;
        }
       // call_pyhepc( 1 );
    }
}