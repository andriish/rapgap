C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     1-LOOP FORM FACTORS, WEAK CORRECTIONS
C
C     FOR ELECTRON QUARK SCATTERING:
C      DIRECT BOX:  (T,S)  IN CFWG1
C     CROSSED BOX:  (T,U)  IN CFWG2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSWG1L(GS,GT,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSWG1L,CM2,HSCIR,HSCWLL,HSCWQL
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      HSWG1L=GT/2D0*HSCIR(GT,MU2)
     *       +CM2/2D0*(HSCWLL(GS,CM2)+HSCWQL(GS,CM2,MU2))
      HSWG1L=HSWG1L*ALP2PI
      END
