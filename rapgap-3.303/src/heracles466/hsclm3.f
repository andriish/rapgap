C=======================================================================
      FUNCTION HSCLM3(Q2,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSCLM3,CM2,W,HSCLN,CHI,CLH
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      W=CM2/DCMPLX(Q2,0D0)
      IF(REAL(W).GT.0) GOTO 100
      CHI=SQRT(1.-4.*W)
      CLH=HSCLN((1+CHI)/(CHI-1))
      HSCLM3=5./6.-2./3.*W+(2.*W+1.)/3.*CHI*CLH
     1        +2./3.*W*(W+2.)*CLH*CLH
      GOTO 400
100   IF (Q2 - 4.*REAL(CM2)) 200,300,300
200   CHI=SQRT(4.*W-1)
      CLH=ATAN(1./DREAL(CHI))
      HSCLM3=5./6.-2./3.*W+2./3.*(2.*W+1)*CHI*CLH
     1        -8./3.*W*(W+2.)*CLH*CLH
      GOTO 400
300   CHI=SQRT(1.-4.*W)
      CLH=HSCLN( (1.+CHI)/(1.-CHI) )
      HSCLM3=5./6.-2.*W/3.+(2.*W+1.)/3.*CHI*CLH
     1       +2./3.*W*(W+2.)*(CLH*CLH - PI*PI)
     2       -DCMPLX(0D0,1D0)*PI*( (2.*W+1.)/3.*CHI+2./3.*W*(W+2.)*CLH)
400   CONTINUE
      RETURN
      END
