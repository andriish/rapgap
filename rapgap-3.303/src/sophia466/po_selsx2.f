

      SUBROUTINE PO_SELSX2(XS1,XS2,XMIN,XMAX,AS1,AS2,IREJ)
C***********************************************************************
C
C    select x values of soft string ends using PO_RNDBET
C
C    (taken from PHOJET v1.12, R.E. 08/98)
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      SAVE

      DIMENSION XS1(2),XS2(2)
      DIMENSION XMIN(2),XMAX(2)

      IREJ = 0

      GAM1 = +1.5D0 + 1.D0
      GAM2 = -0.5D0 + 1.D0
      BET1 = -0.5D0 + 1.D0
      BET2 = -0.5D0 + 1.D0

      ITRY0 = 0
      DO 100 I=1,100

        ITRY1 = 0
 10     CONTINUE
          X1 = PO_RNDBET(GAM1,BET1)
          ITRY1 = ITRY1+1
          IF(ITRY1.GE.50) THEN
            IREJ = 1
            RETURN
          ENDIF
        IF((X1.LE.XMIN(1)).OR.(X1.GE.XMAX(1))) GOTO 10

        ITRY2 = 0
 11     CONTINUE
          X2 = PO_RNDBET(GAM2,BET2)
          ITRY2 = ITRY2+1
          IF(ITRY2.GE.50) THEN
            IREJ = 2
            RETURN
          ENDIF
        IF((X2.LE.XMIN(2)).OR.(X2.GE.XMAX(2))) GOTO 11

        X3 = 1.D0 - X1
        X4 = 1.D0 - X2
        IF(X1*X2.GT.AS1) THEN
          IF(X3*X4.GT.AS2) GOTO 300
        ENDIF
        ITRY0 = ITRY0+1

 100  CONTINUE

      IREJ = 3
      RETURN

 300  CONTINUE

      XS1(1) = X1
      XS1(2) = X3

      XS2(1) = X2
      XS2(2) = X4

      END
