


      SUBROUTINE PO_MSHELL(PA1,PA2,XM1,XM2,P1,P2)
C********************************************************************
C
C    rescaling of momenta of two partons to put both
C                                       on mass shell
C
C    input:       PA1,PA2   input momentum vectors
C                 XM1,2     desired masses of particles afterwards
C                 P1,P2     changed momentum vectors
C
C    (taken from PHOJET 1.12, R.E. 08/98)
C
C********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      PARAMETER ( DEPS = 1.D-5 )

      DIMENSION PA1(4),PA2(4),P1(4),P2(4)

C  Lorentz transformation into system CMS
      PX = PA1(1)+PA2(1)
      PY = PA1(2)+PA2(2)
      PZ = PA1(3)+PA2(3)
      EE = PA1(4)+PA2(4)
      XMS = EE**2-PX**2-PY**2-PZ**2
      XMS = SQRT(XMS)
      BGX = PX/XMS
      BGY = PY/XMS
      BGZ = PZ/XMS
      GAM = EE/XMS
      CALL PO_ALTRA(GAM,-BGX,-BGY,-BGZ,PA1(1),PA1(2),PA1(3),
     &           PA1(4),PTOT1,P1(1),P1(2),P1(3),P1(4))
C  rotation angles
      PTOT1 = MAX(DEPS,PTOT1)
      COD= P1(3)/PTOT1
      SID= SQRT((1.D0-COD)*(1.D0+COD))
      COF=1.D0
      SIF=0.D0
      IF(PTOT1*SID.GT.1.D-5) THEN
        COF=P1(1)/(SID*PTOT1)
        SIF=P1(2)/(SID*PTOT1)
        ANORF=SQRT(COF*COF+SIF*SIF)
        COF=COF/ANORF
        SIF=SIF/ANORF
      ENDIF

C  new CM momentum and energies (for masses XM1,XM2)
      XM12 = XM1**2
      XM22 = XM2**2
      SS   = XMS**2
      PCMP = PO_XLAM(SS,XM12,XM22)/(2.D0*XMS)
      EE1  = SQRT(XM12+PCMP**2)
      EE2  = XMS-EE1
C  back rotation
      CALL PO_TRANS(0.D0,0.D0,PCMP,COD,SID,COF,SIF,XX,YY,ZZ)
      CALL PO_ALTRA(GAM,BGX,BGY,BGZ,XX,YY,ZZ,EE1,
     &           PTOT1,P1(1),P1(2),P1(3),P1(4))
      CALL PO_ALTRA(GAM,BGX,BGY,BGZ,-XX,-YY,-ZZ,EE2,
     &           PTOT2,P2(1),P2(2),P2(3),P2(4))

      END
