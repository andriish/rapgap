      function getorig(a)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer a,getorig
      getorig=K(a,3)
      end 

      function getq2q()
#include "rgpara1.inc"
      double precision getq2q
      getq2q=Q2Q
      end 