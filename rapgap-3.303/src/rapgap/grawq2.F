*#**********************************************************************
*#
*#    SUBROUTINE GRAWQ2
*#
*# PURPOSE: Decide whether to write event for weighted files
*#
*# INPUT: RAPGAP Common RGPARAS
*#
*# OUTPUT: IWRITE = 0 If event to be rejected
*#                = 1 If event to be kept
*#
*# CALLED BY: GRAPGA
*#
*# CALLING: draprn()
*#
*# AUTHOR: Julian Phillips                 CREATED AT: 96/10/01
*#**********************************************************************
      SUBROUTINE GRAWQ2(IWRITE,Q2EG)
	IMPLICIT NONE
#include "rgparas.inc"
#include "rgpara1.inc"
      Integer NCHK
      PARAMETER (NCHK=1000)
      INTEGER  IWRITE
      Double Precision draprn
	Real Q2EG
      EXTERNAL draprn
      DOUBLE PRECISION WTEST

      INTEGER ILOOP,ICALL,IKEEP
      DATA ILOOP/0/, ICALL/0/, IKEEP/0/

      IWRITE = 0
      ICALL = ICALL + 1
      ILOOP = ILOOP + 1

* Contruct weight
      W_Q2 = 1.D0
      IF(Q2.GT.0.D0) W_Q2 = DMAX1(1.D0,Q2START/DBLE(Q2EG))

* Keep fraction 1/W_Q2 of sample, weighted by W_Q2
      WTEST = draprn()
      IF(WTEST.LT.1.D0/W_Q2) THEN
         IWRITE=1
         IKEEP = IKEEP + 1
      ENDIF

      IF(ILOOP.EQ.NCHK) THEN
         ILOOP = 0
         WRITE(6,*)'Weighting: from ',ICALL,' kept ',IKEEP
      ENDIF

      RETURN
      END
