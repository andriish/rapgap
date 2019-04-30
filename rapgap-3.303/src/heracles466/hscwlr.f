C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCWLR(FS,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSCWLR,CM2,CSX,HSCLN,HSSPEN
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      CSX=DCMPLX(FS,1D-6)
      HSCWLR=-1D0/FS*(HSCLN(-CSX/CM2)*HSCLN(-CSX/CM2)/2D0
     .                +HSSPEN((CSX-CM2)/CSX)-PI*PI/3D0
     .                +DCMPLX(0D0,1D0)*PI*HSCLN((CSX-CM2)/CSX))
      RETURN
      END
