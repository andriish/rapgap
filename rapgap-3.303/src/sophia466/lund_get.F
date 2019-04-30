

      SUBROUTINE lund_get(I,IFL,PX,PY,PZ,EE,XM)
C***********************************************************************
C
C   read final states from Lund common block
C
C                                                      (R.E. 08/98)
C
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

#include "rglujets.inc"
#include "rgludat1.inc"
#include "rgludat3.inc"

      PX = PYP(I,1)
      PY = PYP(I,2)
      PZ = PYP(I,3)
      EE = PYP(I,4)
      XM = PYP(I,5)

      Il = PYK(I,8)

C  convert particle ID

      IFL = ICON_PDG_SIB(Il)

      END
