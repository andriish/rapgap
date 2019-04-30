

      SUBROUTINE DECPAR(LA,P0,ND,LL,P)
C***********************************************************************
C
C   This subroutine generates the decay of a particle
C   with ID = LA, and 5-momentum P0(1:5)
C   into ND particles of 5-momenta P(j,1:5) (j=1:ND)
C 
C   If the initial particle code is LA=0
C   then ND and LL(1:ND) are considered as  input and
C   the routine generates a phase space decay into ND
C   particles of codes LL(1:nd)
C
C   (taken from SIBYLL 1.7, muon decay corrected, R.E. 04/98)
C 
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       COMMON /S_CSYDEC/ CBR(102), IDB(49), KDEC(612), LBARP(49)
      COMMON /S_MASS1/ AM(49), AM2(49)
      SAVE

      DIMENSION P0(5), LL(10), P(10,5)
      DIMENSION PV(10,5), RORD(10), UE(3),BE(3), FACN(3:10)

      DATA FACN /2.D0,5.D0,15.D0,60.D0,250.D0,
     +          1500.D0,12000.D0,120000.D0/
      DATA PI /3.1415926D0/

C...c.m.s. Momentum in two particle decays
      PAWT(A,B,C) = SQRT((A**2-(B+C)**2)*(A**2-(B-C)**2))/(2.D0*A)

C...Phase space decay into the particles in the list
      IF (LA .EQ. 0)  THEN
          MAT = 0
          MBST = 0
          PS = 0.
          DO J=1,ND
             P (J,5) = AM(IABS(LL(J)))
             PV(J,5) = AM(IABS(LL(J)))
             PS = PS+P(J,5)
          ENDDO
          DO J=1,4
             PV(1,J) = P0(J)
          ENDDO
          PV(1,5) = P0(5)
          GOTO 140
      ENDIF

C...Choose decay channel
      L = IABS(LA)
      ND=0
      IDC = IDB(L)-1
      IF (IDC+1 .LE.0)  RETURN
      RBR = RNDM(0)
110   IDC=IDC+1
      IF(RBR.GT.CBR(IDC))  GOTO 110

      KD =6*(IDC-1)+1
      ND = KDEC(KD)
      MAT= KDEC(KD+1)

      MBST=0
      IF (MAT .GT.0 .AND. P0(4) .GT. 20.D0*P0(5)) MBST=1
      IF (MAT .GT.0 .AND. MBST .EQ. 0)
     +        BETA = SQRT(P0(1)**2+P0(2)**2+P0(3)**2)/P0(4)

      PS = 0.D0
      DO J=1,ND
         LL(J) = KDEC(KD+1+J)
         P(J,5)  = AM(LL(J))
         PV(J,5) = AM(LL(J))
         PS = PS + P(J,5)
      ENDDO
      DO J=1,4
         PV(1,J) = 0.D0
         IF (MBST .EQ. 0)  PV(1,J) = P0(J)
      ENDDO
      IF (MBST .EQ. 1)  PV(1,4) = P0(5)
      PV(1,5) = P0(5)

140   IF (ND .EQ. 2) GOTO 280

      IF (ND .EQ. 1)  THEN
         DO J=1,4
            P(1,J) = P0(J)
         ENDDO
         RETURN
      ENDIF

C...Calculate maximum weight for ND-particle decay
      WWTMAX = 1.D0/FACN(ND)
      PMAX=PV(1,5)-PS+P(ND,5)
      PMIN=0.D0
      DO IL=ND-1,1,-1
         PMAX = PMAX+P(IL,5)
         PMIN = PMIN+P(IL+1,5)
         WWTMAX = WWTMAX*PAWT(PMAX,PMIN,P(IL,5))
      ENDDO

C...generation of the masses, compute weight, if rejected try again
240   RORD(1) = 1.D0
      DO 260 IL1=2,ND-1
      RSAV = RNDM(0)
      DO 250 IL2=IL1-1,1,-1
      IF(RSAV.LE.RORD(IL2))   GOTO 260
250     RORD(IL2+1)=RORD(IL2)
260     RORD(IL2+1)=RSAV
      RORD(ND) = 0.D0
      WT = 1.D0
      DO 270 IL=ND-1,1,-1
      PV(IL,5)=PV(IL+1,5)+P(IL,5)+(RORD(IL)-RORD(IL+1))*(PV(1,5)-PS)
270   WT=WT*PAWT(PV(IL,5),PV(IL+1,5),P(IL,5))
      IF (WT.LT.RNDM(0)*WWTMAX)   GOTO 240

C...Perform two particle decays in respective cm frame
280   DO 300 IL=1,ND-1
      PA=PAWT(PV(IL,5),PV(IL+1,5),P(IL,5))
      UE(3)=2.D0*RNDM(0)-1.D0
      PHI=2.D0*PI*RNDM(0)
      UT = SQRT(1.D0-UE(3)**2)
      UE(1) = UT*COS(PHI)
      UE(2) = UT*SIN(PHI)
      DO 290 J=1,3
      P(IL,J)=PA*UE(J)
290   PV(IL+1,J)=-PA*UE(J)
      P(IL,4)=SQRT(PA**2+P(IL,5)**2)
300   PV(IL+1,4)=SQRT(PA**2+PV(IL+1,5)**2)

C...Lorentz transform decay products to lab frame
      DO 310 J=1,4
310   P(ND,J)=PV(ND,J)
      DO 340 IL=ND-1,1,-1
      DO 320 J=1,3
320   BE(J)=PV(IL,J)/PV(IL,4)
      GA=PV(IL,4)/PV(IL,5)
      DO 340 I=IL,ND
      BEP = BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
      DO 330 J=1,3
330   P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J)
340   P(I,4)=GA*(P(I,4)+BEP)

C...Weak decays
        IF (MAT .EQ. 1)  THEN
           F1=P(2,4)*P(3,4)-P(2,1)*P(3,1)-P(2,2)*P(3,2)-P(2,3)*P(3,3)   
           IF (MBST.EQ.1)  WT = P0(5)*P(1,4)*F1
           IF (MBST.EQ.0)  
     +     WT=F1*(P(1,4)*P0(4)-P(1,1)*P0(1)-P(1,2)*P0(2)-P(1,3)*P0(3))
           WTMAX = P0(5)**4/16.D0
           IF(WT.LT.RNDM(0)*WTMAX)   GOTO 240
        ENDIF


C...Boost back for rapidly moving particle
      IF (MBST .EQ. 1)   THEN
         DO 440 J=1,3
440      BE(J)=P0(J)/P0(4)
         GA= P0(4)/P0(5)
         DO 460 I=1,ND
         BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
         DO 450 J=1,3
450         P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J)
460         P(I,4)=GA*(P(I,4)+BEP)
      ENDIF

C...labels for antiparticle decay
      IF (LA .LT. 0 .AND. L .GT. 18)  THEN
           DO J=1,ND
            LL(J) = LBARP(LL(J))
         ENDDO
      ENDIF

      END
