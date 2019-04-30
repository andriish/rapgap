

      SUBROUTINE lund_put(I,IFL,PX,PY,PZ,EE)
C***********************************************************************
C
C   store initial configuration into Lund common block
C
C                                                      (R.E. 08/98)
C
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

#include "rglujets.inc"
#include "rgludat1.inc"
#include "rgludat3.inc"

      if(IFL.eq.1) then
        Il = 2
      else if(IFL.eq.-1) then
        Il = -2
      else if(IFL.eq.2) then
        Il = 1
      else if(IFL.eq.-2) then
        Il = -1
      else if(IFL.eq.11) then
        Il = 2203
      else if(IFL.eq.12) then
        Il = 2101
      else if(IFL.eq.21) then
        Il = 2103
      else if(IFL.eq.22) then
        Il = 1103
      else
        print *,' lund_put: unkown flavor code',IFL
      endif

      P(I,1) = PX
      P(I,2) = PY
      P(I,3) = PZ
      P(I,4) = EE
      P(I,5) = (EE-PZ)*(EE+PZ)-PX**2-PY**2

      K(I,1) = 1
      K(I,2) = Il
      K(I,3) = 0
      K(I,4) = 0
      K(I,5) = 0

      N = I

      END
