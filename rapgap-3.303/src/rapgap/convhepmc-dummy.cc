//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <cmath>
#include <iostream>

using namespace std;

extern "C" {   
int ncount ;
    void convhepmc_(int & ievent, int & iproc, double & xsec, double & xsece,
                    int& flav1, int& flav2,double &  x1,double &  x2,double &  q2pdfeval,
                    double &  xf1mom, double & xf2mom ,int&  pdf1,int& pdf2
                   )
                   {
      if ( ncount < 10) {
          cout << " CASCADE: dummy version of convhepmc is used "<< endl;    
          ++ncount;
      }
   }   
}
