C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCIR(FS,MQ2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSCIR
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      IF (FS) 1,2,3
    1 DLNE=DLOG(-ME2/FS)
      DLNQ=DLOG(-MQ2/FS)
      AIR=-1D0/FS*(DLNE*DLNE/4D0+DLNQ*DLNQ/4D0+PI*PI/6D0)
      HSCIR=DCMPLX(AIR,0D0)
      RETURN
    2 HSCIR=(0D0,0D0)
      RETURN
    3 DLNE=DLOG(ME2/FS)
      DLNQ=DLOG(MQ2/FS)
      AIR=-1D0/FS*(DLNE*DLNE/4D0+DLNQ*DLNQ/4D0+2D0*PI*PI/3D0)
      HSCIR=DCMPLX(AIR,0D0)
      RETURN
      END
