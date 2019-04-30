C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE FIUSER(X,Q2,ZF1,ZF2,IPDFR)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      ZF1=0D0
      ZF2=0D0
      IPDFR=0

C--->
C    insert here your conditions on x,Q2 for range of applicability 
C    of your parametrization. If outside, set IPDFR = 1 and return. 
C    The calling routine will then use parton distributions as specified
C    by ILIB and ICODE to calculate structure functions
C 
c      IF (X.LT.XFMIN.OR.X.GT.XFMAX.
c     &.OR.Q2.LT.QFMIN.OR.Q2.GT.QFMAX) THEN
c        IPDFR=1
c        RETURN
c      ENDIF 
C--->

C--->
C    insert here your parametrization for F1 and F2 (electromagnetic
C    part, no contribution from Z exchange).
C    Note: F1 = (F2-FL)/(2*X)
C
c      ZF1=
c      ZF2=
C--->
      RETURN
      END
