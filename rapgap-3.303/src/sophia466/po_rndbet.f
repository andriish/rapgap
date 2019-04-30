

      DOUBLE PRECISION FUNCTION PO_RNDBET(GAM,ETA)
C********************************************************************
C
C    random number generation from beta
C    distribution in region  0 < X < 1.
C    F(X) = X**(GAM-1.)*(1.-X)**(ETA-1)*GAMM(ETA+GAM) / (GAMM(GAM
C                                                        *GAMM(ETA))
C
C    (taken from PHOJET v1.12, R.E. 08/98)
C
C********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      SAVE

      Y = PO_RNDGAM(1.D0,GAM)
      Z = PO_RNDGAM(1.D0,ETA)
      PO_RNDBET = Y/(Y+Z)

      END
