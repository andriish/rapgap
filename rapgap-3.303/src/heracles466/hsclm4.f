C=======================================================================
      FUNCTION HSCLM4(Q2,CM1,CM2,MF2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSCLM4,CM1,CM2
     1          ,W1,W2,X1,X2,C12,HSCLN,HSSPEN
      RM1=SQRT(DREAL(CM1))
      RM2=SQRT(DREAL(CM2))
      SM2=(RM1+RM2)*(RM1+RM2)
      DM2=(RM1-RM2)*(RM1-RM2)
      W1=CM1/DCMPLX(Q2,0D0)
      W2=CM2/DCMPLX(Q2,0D0)
      IF (REAL(CM1).EQ.0) GOTO 401
      IF (REAL(CM2).EQ.0) GOTO 402
      IF ( (Q2.GT.DM2).AND.(Q2.LT.SM2) ) GOTO 100
      X1=(1.-W1+W2)/2.+SQRT((1.-W1+W2)*(1.-W1+W2)-4.*W2)/2.
      X2=(1.-W1+W2)/2.-SQRT((1.-W1+W2)*(1.-W1+W2)-4.*W2)/2.
      GOTO 200
100   X1=(1.-W1+W2)/2.+DCMPLX(0D0,1D0)/2.
     1                 *SQRT(4.*W2-(1.-W1+W2)*(1.-W1+W2) )
      X2=(1.-W1+W2)/2.-DCMPLX(0D0,1D0)/2.
     1                 *SQRT(4.*W2-(1.-W1+W2)*(1.-W1+W2) )
200   C12=HSCLN(CM1/CM2)/2.
      HSCLM4=1./6. + (W1+W2)/(W1-W2)*C12 - (W1-W2)/3.*C12
     1       +(W1+W2+1.)/3.*(C12-1.)
     2       +(W1+W2+1.)/3.*(X1*HSCLN(X1/(X1-1.))
     3                       + X2*HSCLN(-X2/(1.-X2)) )
     4       -2./3.*(W1+W2+W1*W2)*HSCLN(X1/(X1-1.))*HSCLN(-X2/(1.-X2))
      GOTO 500
401   W1=W2
      X1=HSCLN(CM2/MF2)
      X2=HSCLN((CM2-Q2)/MF2)
      GOTO 404
402   X1=HSCLN(CM1/MF2)
      X2=HSCLN((CM1-Q2)/MF2)
404   HSCLM4=.5 - W1/3. + (1.-W1*W1)/3.*HSCLN((W1-1.)/W1)
     1       +2./3.*X1 + W1/3.*(X2*X2-X1*X1+2.*HSSPEN(1./(1.-W1)) )
500   RETURN
      END
