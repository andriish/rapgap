C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C---FILL TERMINATING RECORD OF COMMON HEPEVT
C
      SUBROUTINE HSTERM
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C---------------------------------------------------------------------
      COMMON /HSNUME/ SIGTOT,SIGTRR,SIGG(20),SIGGRR(20),NEVENT,NEVE(20)
      COMMON /HSGSW/  SW,CW,SW2,CW2
     +              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     +              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /PYSTFUC/ PYSTOP,PYSLAM,NPYMOD,NPYMAX,NPYMIN
      REAL           PYSTOP,PYSLAM
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
#include "hepevt.inc"
      NEVHEP=-2

C---CROSS SECTIONS
      PHEP(1,1)=SIGTOT/1D3
      PHEP(1,2)=SIGTRR/1D3
      DO 10 I=1,20
         PHEP(1,2+I)=SIGG(I)/1D3
   10 PHEP(1,22+I)=SIGGRR(I)/1D3
C---MASSES, SW2
      PHEP(1,43)=SW2
      PHEP(1,44)=MW
C---EVENT NUMBERS
      ISTHEP(1)=NEVENT
      DO 20 I=1,20
   20 ISTHEP(1+I)=NEVE(I)
      NHEP=44

      RETURN
      END
