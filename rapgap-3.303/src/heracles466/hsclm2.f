C=======================================================================
      FUNCTION HSCLM2(Q2,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSCLM2,W,HSCLN,HSSPEN,CM2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      W=CM2/DCMPLX(Q2,0D0)
      IF(REAL(W).GT.0) GOTO 100
      HSCLM2=.5-2.*(2.+W)-(2.*W+3)*HSCLN(-W)+
     1        2.*(1+W)**2*(HSSPEN(1+1/W)-PI*PI/6)
      GOTO 200
100   HSCLM2=.5-2.*(2.+W)-(2.*W+3)*HSCLN(W)+
     1        2.*(1+W)**2*(HSCLN(W)*HSCLN((W+1.)/W)-HSSPEN(-1/W))
     2        -DCMPLX(0D0,1D0)*PI*(3.+2.*W-2.*(1+W)**2*HSCLN((1.+W)/W))
200   CONTINUE
      RETURN
      END
