C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE BRASSE(ZT,ZAMF2,ZF1,ZF2)
      IMPLICIT DOUBLE PRECISION (Z)
      COMMON /HSSINC/ S,AP,AMP2,W2PIT,W2MIN,W2TR
C
      FMF2=SNGL(ZAMF2)
      IF(FMF2.LE.W2MIN)   CALL WABC(W2MIN)
      IF(FMF2.GT.W2MIN.AND.FMF2.LE.W2TR) CALL WABC(FMF2)
      IF(FMF2.GT.W2TR )   CALL WABC(W2TR )
      T=SNGL(ZT)
      CALL RF12(T,FMF2,RF1,RF2)
      ZF1=DBLE(RF1)
      ZF2=DBLE(RF2)
      END
