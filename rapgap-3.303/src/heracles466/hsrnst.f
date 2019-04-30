C===========================================================
      SUBROUTINE HSRNST(NA1,NA2,NA3,NB1)
C     REAL*4 U,C,CD,CM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /HSRNDC/ U(97),C,CD,CM,I,J
C
      MA1 = NA1
      MA2 = NA2
      MA3 = NA3
      MB1 = NB1
      I   = 97
      J   = 33
      DO 20 II2 = 1,97
        S = 0D0
        T = 0.5D0
        DO 10 II1 = 1,24
          MAT  = MOD(MOD(MA1*MA2,179)*MA3,179)
          MA1  = MA2
          MA2  = MA3
          MA3  = MAT
          MB1  = MOD(53*MB1+1,169)
          IF ( MOD(MB1*MAT,64).GE.32 ) S = S+T
10        T    = 0.5D0*T
20      U(II2) = S
      C  =   362436.D0/16777216.D0
      CD =  7654321.D0/16777216.D0
      CM = 16777213.D0/16777216.D0
      RETURN
      END
