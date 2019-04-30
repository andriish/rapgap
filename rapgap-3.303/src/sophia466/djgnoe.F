
C***********************************************************************

      SUBROUTINE DJGNOE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /SPPASS/ NSOPH,NSPOUT,NFAILP,NSPACC
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      REAL*4 CUT,PARL,X,Y,W2,Q2,U
      SAVE
#include "rgluco.inc"
      N=0
      LST(21)=99
      NSPOUT=NSPOUT+1
      RETURN
      END
