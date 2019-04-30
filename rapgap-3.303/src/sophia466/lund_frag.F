

      SUBROUTINE lund_frag(SQS,NP)
C***********************************************************************
C
C   interface to Lund/Jetset fragmentation
C
C                                    (R.E. 08/98)
C
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

#include "rglujets.inc"
#include "rgludat1.inc"
#include "rgludat3.inc"


C  no title page

c        MSTU(12) = 0

C  define some particles as stable

        MSTJ_22 = MSTJ(22)
        MSTJ(22) = 2

C  in addition pi0 stable

        KC=PYCOMP(111)
        MDCY_pi0 = MDCY(KC,1)
        MDCY(KC,1)=0


C  switch popcorn effect off

        MSTJ_12 = MSTJ(12)
        MSTJ(12) = 1

C  suppress all warning and error messages

        MSTU_22 = MSTU(22)
        MSTU(22) = 0
        MSTU_25 = MSTU(25)
        MSTU(25) = 0


C  set energy dependent parameters

      PARJ_36 = PARJ(36)
      IF(SQS.LT.2.D0) THEN
        PARJ(36) = 0.1D0
      ELSE IF(SQS.LT.4.D0) THEN
        PARJ(36) = 0.7D0*(SQS-2.D0)/2.D0+0.1D0
      ELSE
        PARJ(36) = 0.8D0
      ENDIF

C  fragment string configuration

      II = MSTU(21)
      MSTU(21) = 1
      CALL PYEXEC
      MSTU(21) = II
      
C  event accepted?

      if(MSTU(24).ne.0) then
        NP = -1
        return
      endif

      CALL PYEDIT(1)

      NP = PYP(0,1)

C  restore original parameter settings

      MSTJ(12) = MSTJ_12
      MSTJ(22) = MSTJ_22

      MSTU(22) = MSTU_22
      MSTU(25) = MSTU_25

      MDCY(KC,1) = MDCY_pi0

      PARJ(36) = PARJ_36

      END
