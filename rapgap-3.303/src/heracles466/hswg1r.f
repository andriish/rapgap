C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSWG1R(GS,GT,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSWG1R,CM2,HSD13C,HSCMW,HSCWLR,HSCWQR
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      GU=-(GS+GT)
      HSWG1R=-GU*HSD13C(GT,GS,CM2)+GT/2D0*HSCMW(GT,CM2)
     *       +CM2/2D0*(HSCWLR(GS,CM2)+HSCWQR(GS,CM2,MU2))
      HSWG1R=HSWG1R*ALP2PI
      END
