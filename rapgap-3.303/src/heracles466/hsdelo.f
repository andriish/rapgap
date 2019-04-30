
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE HSDELO(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      IF (IOPEGM.GT.0) THEN
        DELTA=EGMIN
        RETURN
      ELSE
        SIGMA=EPRO/EELE
        XMY=Y-(1D0-X*Y)*SIGMA
        XPY=Y+(1D0-X*Y)*SIGMA
        OMEGA=2D0*EELE*SIGMA*Y*(1D0-X)
     *          /(XPY+DSQRT(4D0*X*Y*SIGMA*(1D0-Y)+XMY*XMY))
        EQUA=X*SIGMA*EELE
        EES=EELE*(1D0-Y+X*Y*SIGMA)
        DELTA=DMIN1(EELE,EES,OMEGA)*DELEPS
      ENDIF
      END
