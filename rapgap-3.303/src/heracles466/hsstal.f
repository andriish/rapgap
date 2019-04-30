
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   Interface to the parametrization of ALLM for F2 (1997 update)
C   H.Abramowicz, A.Levy, DESY 97-251 (Dec.1997) (hep-ph/9712415)
C
      SUBROUTINE HSSTAL(X,Q2,F1R,F2R)
      DOUBLE PRECISION X,Q2,F1R,F2R
      SX=SNGL(X)
      SQ2=SNGL(Q2)
      F2E=F2ALLM(SX,SQ2)
      F2R=DBLE(F2E)
      F1R=F2R/2D0/X
      RETURN
      END
