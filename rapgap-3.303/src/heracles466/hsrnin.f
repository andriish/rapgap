C==========================================================
      SUBROUTINE HSRNIN(UIN,CIN,CDIN,CMIN,IIN,JIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     REAL*4 UIN(97),CIN,CDIN,CMIN
C     REAL*4 U,C,CD,CM
C
      COMMON /HSRNDC/ U(97),C,CD,CM,I,J
      DIMENSION UIN(97)
C
      DO 10 KKK = 1,97
10      U(KKK) = UIN(KKK)
      C  = CIN
      CD = CDIN
      CM = CMIN
      I  = IIN
      J  = JIN
      RETURN
      END
