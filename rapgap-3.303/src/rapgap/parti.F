      SUBROUTINE PARTI(IP1,YX,YVAL,WEIGHT,IFL,IST)
      IMPLICIT NONE
	Integer IP1,IFL,IST
	Double Precision YX,YVAL,WEIGHT
#include "rglujets.inc"
#include "rgluco.inc"
#include "rgpypars.inc"
#include "rgpara1.inc"
#include "rgpara.inc"
#include "rgpart.inc"
#include "rgrapgki.inc"
#include "rgdiffr.inc"
#include "hepevt.inc"
#include "rghs45.inc"
C     from HERACLES;
	Integer ICHNN,LST,IRES
      COMMON /HSCHNN/ ICHNN
      COMMON/EPPARA/LST(30),IRES(2)
	Double Precision PHI,DBGAM
      COMMON/DIFFA/ PHI
      DOUBLE PRECISION ME
      DIMENSION DBGAM(4)
	Double Precision XG1,ALPH_EM,FGAM,SQRTS,EC,X11,COST2,SPHE,CPHE
	Double Precision PEP,PEZ,PDLE,CSPHI,XP1,BOCHCK,THE,PHIE
	Double Precision SPHI,STHETA,PT
	Integer I,NIPH,IE,NFI,IG
C IST = 0, ONLY PARTICLE MOMENTA, BUT SCALE Q2Q IN STRUCTURE FUNCTION
C NOT YET DEFINED.
C IST = 1, CALCULATE ONLY STRUCTURE FUNCTION SINCE Q2Q NOW DEFINED
C
C IF PARTICLE IP1 = 11 (ELECTRON) USE EQUIVALENT PHOTON APPROXIMATION
C
c      write(6,*) 'parti  Q2 = ',Q2,' yx = ',yx
      XG1 = DBLE(XEL)
      YVAL =  999999.D0
      WEIGHT = -999999.D0
      ME = PYMASS(IP1)
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)

C IF 1ST INCOMING PARTICLE IS A ELECTRON AND LOOKING FOR DIRECT PHOTON
      IF((IABS(IP1).EQ.11.OR.IABS(IP1).EQ.13)) THEN
C
         IF(IABS(IP1).EQ.11.OR.IABS(IP1).EQ.13) THEN
C...          WEIZAECKER WILLIAMS APPROXIMATION FOR BREMSTRAHLUNG (EPA)
C...          WITH Q2 DEPENDENCE
            FGAM = (1.D0 - YX +YX*YX/2.D0)/Q2/YX - ME*ME*YX/Q2/Q2

            YVAL = ALPH_EM * FGAM / PI
         ENDIF
         IF(QMI.EQ.0.D0) THEN
            WEIGHT = YX*Q2*DLOG(YMAX/YMIN)*DLOG(Q2MAX/Q2MIN)
         ELSE
c new try weighting with 1/q**4
            WEIGHT = YX*DLOG(YMAX/YMIN)*Q2**2*(Q2MAX-Q2MIN)/Q2MIN/
     +      Q2MAX
         ENDIF
C...

C...      CALCULATE COS OF ELECTRON VIA FORMULA GIVEN IN CAMPBELL ET AL
C...      CAN.J.PHYS. 39 1981 (1742)
C...      VIA X Q2 Y
C...      X11= Q2/(2P.Q)
         SQRTS=DSQRT(SSS)
         EC= 0.5D0*SQRTS
         IF(IHERAC.EQ.0) THEN
            X11= Q2/YX/SSS
            COST2= ((1.D0 - YX) - X11*YX)/((1.D0 - YX) + X11*YX)
            IF(COST2.GT.1.) WRITE(6,*) 'YX Q2 COST2 ',YX,Q2,COST2
C...      COST2 IS SCATTERING ANGLE OF ELECTRON
            THE= DACOS(COST2)
            PHIE = PHI
            IF(PHIE.GT.7.D0) write(6,*) 'fatal PHIE = ',PHIE
            SPHE=DSIN(PHIE)
            CPHE=DCOS(PHIE)
C...      PEL1 = 4 VECTOR OF SCATTERED ELECTRON IN EP CMS SYSTEM
            N=N+1
            K(N,1)=1
            IF(INTER.LT.2) THEN
               K(N,2)=K(1,2)
            ELSEIF(INTER.EQ.2) THEN
               K(N,2)=ISIGN(1,K(1,2))*12
            ENDIF
            K(N,3)=IFL
c. 2*P.l_e
c new gymnastics for numerical stability
            PEP = DSQRT(P(2,3)**2 + P(2,5)**2)
            PEZ = P(2,3)
            PDLE = 2.D0*DOT1(1,2)
            X11 = Q2/PDLE/YX
            P(N,5)=PYMASS(K(N,2))
            P(N,4) = P(1,4) + (PEZ*(Q2+P(N,5)**2)/ (2.D0*P(1,4)) -
     +      YX*PDLE/2.D0)/(PEP+PEZ)
            P(N,3) = P(1,4) - (PEP*(Q2+P(N,5)**2)/ (2.D0*P(1,4)) + YX*
     +      PDLE/2.0D0)/(PEP+PEZ)
            P(N,3)=-P(N,3)
            PT = Q2 - (Q2 + P(N,5)**2)/P(1,4)/(PEP+PEZ)* (YX*PDLE/2.D0 +
     +      P(2,5)**2*(Q2+P(N,5)**2)/4.D0/ P(1,4)/(PEP+PEZ))
            P(N,1)=DSQRT(DMAX1(0.D0,PT))*CPHE
            P(N,2)=DSQRT(DMAX1(0.D0,PT))*SPHE

C...      PPH = 4 VECTOR OF (VIRTUAL) PHOTON IN EP CMS SYSTEM
            N=N+1
            K(N,1)=21
            K(N,2)=KEB
            K(N,3)=IFL
            P(N,1)= P(1,1) - P(3,1)
            P(N,2)= P(1,2) - P(3,2)
            P(N,3)= P(1,3) - P(3,3)
            P(N,4)= P(1,4) - P(3,4)
            P(N,5)= -SQRT(ABS(DOT1(N,N)))
c            write(6,*) 'parti: ',P(N,5)**2,q2
            NIA1 = N
            NIPH = N
         ELSEIF(IHERAC.EQ.1) THEN
	   
c	      write(6,*) ' herac test ', (ISTHEP(i),i=1,5),(IDHEP(i),i=1,5)
            WEIGHT = 1.D0
            IE=N+1
            NFI = N + 2
            IF(PHEP(4,3).NE.0.0) THEN
               IG=N+2
               NFI = N + 3
c radiative gamma
               K(IG,1)=1
               K(IG,2)=KPH
c radiated gamma has origin in e or e'
c               write(6,*) 'parti: ICHNN=',ICHNN
               IF(ICHNN.EQ.6.OR.ICHNN.EQ.12) THEN
                  K(IG,3)=1
               ELSE
                  K(IG,3)=4
               ENDIF
            ENDIF
c scattered electron
            K(IE,1)=1
cc            K(IE,2)=K(1,2)
		K(IE,2)=IDHEP(1)
		IF(INTER.EQ.2) THEN
               K(IE,2)=IDHEP(1)
		ENDIF
            K(IE,3)=IFL
c exchanged  boson
            IF(IABS(IDHEP(1)).EQ.11) THEN
               KEB=KPH
            ELSEIF(IABS(IDHEP(1)).EQ.12) THEN
               KEB=-24*ISIGN(1,K(1,2))
            ELSE
               write(6,*) ' PARTI: scattered lepton not known ',IDHEP(1)
            ENDIF
            K(NFI,1)=21
            K(NFI,2)=KEB
            K(NFI,3)=IFL
            DO 10  I=1,5
               P(IE,I) = PHEP(I,1)
               P(NFI,I) = PHEP(I,5) - PHEP(I,1) - PHEP(I,3)
               IF(PHEP(4,3).NE.0.0) THEN
                  P(IG,I) = PHEP(I,3)
               ENDIF
   10       CONTINUE
c calculate phi
            CSPHI = P(IE,1)/SQRT(P(IE,1)**2+P(IE,2)**2)
            PHI = DACOS(CSPHI)
            N = NFI
            NIA1 = NFI
            NIPH = NIA1
            DO 20  I = 3,N
   20       P(I,3) = -P(I,3)
            P(NFI,5)= -SQRT(ABS(DOT1(NFI,NFI)))
c         write(6,*) ' parti: q2hs =',q2hs,p(nfi,5)**2
C BOOST TO EP CMS
            CALL PYROBO(3,N,0.D0,0.D0,-CM(1)/CM(4),-CM(2)/CM(4),-CM(3)/
     +      CM(4))

         ELSE
            write(6,*) ' not valid for IHERAC ',IHERAC
         ENDIF
         IF(IRES(1).EQ.0) THEN
            NIR1 = -99999
         ELSEIF(IRES(1).EQ.1) THEN
c          write(6,*) 'parti yx = ',yx,' xel ',xg1
C         XG1 = X_GLUON OF THE PHOTON = E_GLUON/E_ELECTRON
C         YX  = E_PHOTON/E_ELECTRON
C         XP1 = E_GLUON/E_PHOTON
            XP1=XG1/YX
C...      PPH = NOW 4 VECTOR OF HADR COMP OF PHOTON IN EP CMS SYSTEM
            N = N + 1
            K(N+1,1)=21
            K(N+1,2)=21
            K(N+1,3)=IFL
            P(N+1,1)= P(N-1,1)*XP1
            P(N+1,2)= P(N-1,2)*XP1
            P(N+1,3)= P(N-1,3)*XP1
cccc            P(N+1,4)= DABS(P(N-1,3))*XP1
            P(N+1,4)= DABS(P(N-1,4))*XP1
            P(N+1,5)= PYMASS(21)
            NIA1 = N+1
C...      PHA = NOW 4 VECTOR OF PHOTON REMNANT  IN EP CMS SYSTEM
            K(N,1)=1
            K(N,2)=21
            K(N,3)=IFL
            P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
            P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
            P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
            P(N,5)= PYMASS(21)
            P(N,4)= DSQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)
            P(N+1,4)= P(1,4) - P(3,4) - P(N,4)
            NIR1 = N

c            write(6,*) ' parti yx,xg1,xp1',yx,xg1,xp1
c            call dulist(1)
C NOW BOOST TO GAMMA PROTON FRAME
            DBGAM(1) = P(NIPH,1) + P(2,1)
            DBGAM(2) = P(NIPH,2) + P(2,2)
            DBGAM(3) = P(NIPH,3) + P(2,3)
            DBGAM(4) = P(NIPH,4) + P(2,4)
            BOCHCK = (DBGAM(1)/DBGAM(4))**2 + (DBGAM(2)/DBGAM(4))**2 +
     +      (DBGAM(3)/DBGAM(4))**2
            BOCHCK = DSQRT(BOCHCK)
            IF(BOCHCK.GT.0.99999999D0) goto 30
            CALL PYROBO(0,N,0.D0,0.D0,-DBGAM(1)/DBGAM(4),-DBGAM(2)/
     +      DBGAM(4), -DBGAM(3)/DBGAM(4))
            SPHI = PYANGL(P(NIPH,1),P(NIPH,2))
            call PYROBO(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
            STHETA = PYANGL(P(NIPH,3),P(NIPH,1))
            call PYROBO(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
            P(N+1,1) = 0.D0
            P(N+1,2) = 0.D0
            P(N+1,3) = XP1*P(NIPH,3)
            P(N+1,4) = XP1*(P(NIPH,4)*P(2,4)-P(NIPH,3)*P(2,3))/
     +                 (P(2,4)-P(2,3))
            P(N+1,3) = P(N+1,4)
            P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
            P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
            P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
            P(N,4)= P(1,4) - P(N+1,4)- P(3,4)
            P(N,5)= DSQRT(DABS(P(N,1)**2+P(N,2)**2+P(N,3)**2-P(N,4)**2))
            N = N + 1
c            call dulist(1)
            call PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
            call PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
            CALL PYROBO(0,N,0.D0,0.D0,DBGAM(1)/DBGAM(4),DBGAM(2)/
     +      DBGAM(4), DBGAM(3)/DBGAM(4))
c            call dulist(1)
c            pause
         ENDIF
C IF 1ST INCOMING PARTICLE IS A PHOTON
      ELSEIF(IABS(IP1).EQ.22.AND.IRES(1).EQ.0) THEN
C...      4 VECTOR OF (REAL) PHOTON IN CMS SYSTEM
         N=N+1
         K(N,1)=21
         K(N,2)=KPH
         K(N,3)=IFL
         P(N,1)= P(1,1)
         P(N,2)= P(1,2)
         P(N,3)= P(1,3)
         P(N,4)= P(1,4)
         P(N,5)=PYMASS(KPH)
         NIA1 = N
         NIR1 = -99999
         YVAL = 1.D0
         WEIGHT = 1.D0
      ELSEIF(IABS(IP1).EQ.22.AND.IRES(1).EQ.1) THEN
         WEIGHT = 1.0D0
         YVAL = 1.0D0

C         XG1 = X_GLUON OF THE PHOTON = E_GLUON/E_ELECTRON
C         YX  = E_PHOTON/E_ELECTRON
C         XP1 = E_GLUON/E_PHOTON
         XP1=XG1/YX
C...      PPH = NOW 4 VECTOR OF HADR COMP OF PHOTON IN EP CMS SYSTEM
         N = N + 1
         K(N+1,1)=21
         K(N+1,2)=21
         K(N+1,3)=IFL
         P(N+1,1)= P(N-1,1)*XP1
         P(N+1,2)= P(N-1,2)*XP1
         P(N+1,3)= P(N-1,3)*XP1
cccc            P(N+1,4)= DABS(P(N-1,3))*XP1
         P(N+1,4)= DABS(P(N-1,4))*XP1
         P(N+1,5)= PYMASS(21)
         NIA1 = N+1
C...      PHA = NOW 4 VECTOR OF PHOTON REMNANT  IN EP CMS SYSTEM
         K(N,1)=1
         K(N,2)=21
         K(N,3)=IFL
         P(N,1)= P(1,1) - P(N+1,1)- P(3,1)
         P(N,2)= P(1,2) - P(N+1,2)- P(3,2)
         P(N,3)= P(1,3) - P(N+1,3)- P(3,3)
         P(N,5)= PYMASS(21)
         P(N,4)= DSQRT(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)
         P(N+1,4)= P(1,4) - P(3,4) - P(N,4)
         NIR1 = N
         N = N + 1
c            write(6,*) ' parti yx,xg1,xp1',yx,xg1,xp1
c            call dulist(1)
c now look here for protons....
      ELSEIF(IABS(IP1).EQ.2212.AND.IRES(1).EQ.1) THEN
         WEIGHT = YX*DLOG(YMAX/YMIN)
         YVAL = 1.D0/YX
         PHIE = PHI
         IF(PHIE.GT.7.D0) write(6,*) 'fatal PHIE = ',PHIE
         N=N+1
         K(N,1)=1
         K(N,2)=K(1,2)
         K(N,3)=IFL
C...      PPH = 4 VECTOR OF gluon IN PP CMS SYSTEM
         N=N+1
         K(N,1)=21
         K(N,2)=KEB
         K(N,3)=IFL
         P(N,1)= 0
         P(N,2)= 0
         P(N,3)= yx*P(1,3)
         P(N,4)= yx*Abs(P(1,3))
         P(N,5)= -SQRT(ABS(DOT1(N,N)))
c	   write(6,*) KEB,K(N,2)
c            write(6,*) 'parti: ',P(N,5)**2,q2
         NIA1 = N
         P(N-1,1)= P(1,1) - P(NIA1,1)
         P(N-1,2)= P(1,2) - P(NIA1,2)
         P(N-1,3)= P(1,3) - P(NIA1,3)
         P(N-1,4)= P(1,4) - P(NIA1,4)
         P(N-1,5)= -SQRT(ABS(DOT1(N-1,N-1)))
         NIPH = N
ccc         NIR1 = -99999
         NIR1 = 3999

      ENDIF
      RETURN
   30 write(6,*) ' PARTI boost error '
      RETURN
      END
