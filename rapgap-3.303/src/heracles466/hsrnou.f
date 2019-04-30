C==========================================================
      SUBROUTINE HSRNOU(UOUT,COUT,CDOUT,CMOUT,IOUT,JOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     REAL*4 UOUT(97),COUT,CDOUT,CMOUT
C     REAL*4 U,C,CD,CM
      COMMON /HSRNDC/ U(97),C,CD,CM,I,J
      DIMENSION UOUT(97)
C
      DO 10 KKK = 1,97
10      UOUT(KKK) = U(KKK)
      COUT  = C
      CDOUT = CD
      CMOUT = CM
      IOUT  = I
      JOUT  = J
      RETURN
      END
