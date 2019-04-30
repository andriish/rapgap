C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSDELX(XA,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIRCX/ XIRDEL
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      IF (IOPEGM.GT.0) THEN
        DELTA=EGMIN
        RETURN
      ELSE
        X=XA-XIRDEL
        SIGMA=EPRO/EELE
        XMY=Y-(1D0-X*Y)*SIGMA
        XPY=Y+(1D0-X*Y)*SIGMA
        OMEGA=2D0*EELE*SIGMA*Y*(1D0-X)
     *        /(XPY+DSQRT(4D0*X*Y*SIGMA*(1D0-Y)+XMY*XMY))
        EQUA=X*SIGMA*EELE
        EES=EELE*(1D0-Y+X*Y*SIGMA)
C       DELTA=DMIN1(EELE,EES,OMEGA)*DELEPS
C       DELTA=DMIN1(EELE,EES,OMEGA)
        DELTA=OMEGA
      ENDIF
      END
