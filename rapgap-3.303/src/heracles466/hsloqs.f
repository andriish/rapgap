C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSLOQS(Q2,X)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      HSLOQS=1D0-DEXP(-3.37D0*Q2)
      RETURN
      END
