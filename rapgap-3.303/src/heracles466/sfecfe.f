C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SFECFE(SFE,CFE)
      ENTRY         COSI(SFE,CFE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C********************************************************************
C
C     SUBROUTINE OF FLUKA TO GIVE SINE AND COSINE OF AN
C     RANDOM ANGLE UNIFORMLY DISTRIBUTED BETWEEN 0 AND 2*PI
C********************************************************************
C
 10   X=2.0*HSRNDM(-1)-1.0
      Y=HSRNDM(-1)
      X2=X*X
      Y2=Y*Y
      IF (X2+Y2.GT.1.0) GO TO 10
      CFE=(X2-Y2)/(X2+Y2)
      SFE=2.0*X*Y/(X2+Y2)
      RETURN
      END
