

      SUBROUTINE PO_TRANS(XO,YO,ZO,CDE,SDE,CFE,SFE,X,Y,Z)
C**********************************************************************
C
C    rotation of coordinate frame (1) de rotation around y axis
C                                 (2) fe rotation around z axis
C
C    (taken from PHOJET v1.12, R.E. 08/98)
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      X= CDE*CFE*XO-SFE*YO+SDE*CFE*ZO
      Y= CDE*SFE*XO+CFE*YO+SDE*SFE*ZO
      Z=-SDE    *XO       +CDE    *ZO

      END
