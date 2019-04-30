

      DOUBLE PRECISION FUNCTION PO_RNDGAM(ALAM,ETA)
C********************************************************************
C
C    random number selection from gamma distribution
C    F(X) = ALAM**ETA*X**(ETA-1)*EXP(-ALAM*X) / GAM(ETA)
C       
C    (taken from PHOJET v1.12, R.E. 08/98)
C
C********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      SAVE

      NCOU=0
      N = ETA
      F = ETA - N
      IF(F.EQ.0.D0) GOTO 20
   10 R = RNDM(0)
      NCOU=NCOU+1
      IF (NCOU.GE.11) GOTO 20
      IF(R.LT.F/(F+2.71828D0)) GOTO 30
      YYY=LOG(RNDM(0)+1.e-7)/F
      IF(ABS(YYY).GT.50.D0) GOTO 20
      Y = EXP(YYY)
      IF(LOG(RNDM(0)+1.D-7).GT.-Y) GOTO 10
      GOTO 40
   20 Y = 0.D0
      GOTO 50
   30 Y = 1.D0-LOG(RNDM(0)+1.D-7)
      IF(RNDM(0).GT.Y**(F-1.D0)) GOTO 10
   40 IF(N.EQ.0) GOTO 70
   50 Z = 1.D0
      DO 60 I = 1,N
   60 Z = Z*RNDM(0)
      Y = Y-LOG(Z+1.D-7)
   70 PO_RNDGAM = Y/ALAM

      END
