

      DOUBLE PRECISION FUNCTION PO_XLAM(X,Y,Z)
C**********************************************************************
C
C     auxiliary function for two/three particle decay mode
C     (standard LAMBDA**(1/2) function)
C
C     (taken from PHOJET 1.12, R.E. 08/98)
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      YZ=Y-Z
      XLAM=X*X-2.D0*X*(Y+Z)+YZ*YZ
      IF(XLAM.LT.0.D0) XLAM=-XLAM
      PO_XLAM=SQRT(XLAM)

      END
