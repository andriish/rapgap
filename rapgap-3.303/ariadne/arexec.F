C***********************************************************************
C $Id: arexec.F,v 1.1.1.1 2008/11/01 12:35:14 jung Exp $

      SUBROUTINE AREXEC

C...ARiadne subroutine EXECute ariadne

C...The Main driver routine in Ariadne.

      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100,MAXONI=100)
      IMPLICIT DOUBLE PRECISION(A-H, O-P, R-Z)
      IMPLICIT LOGICAL (Q)
      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),QEX(MAXPAR),QQ(MAXPAR),
     $                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),INQ(MAXPAR),
     $                XPMU(MAXPAR),XPA(MAXPAR),PT2GG(MAXPAR),IPART
      SAVE /ARPART/
	REAL PARA
      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/


	COMMON /ARONIA/ IPONI(MAXONI),MEONI(MAXONI),IFLONI(MAXONI),
     $     Q2GONI(MAXONI),PONI(5,MAXONI),NONI,IHQI(MAXSTR,4),NHQ
      SAVE /ARONIA/

      COMMON /PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      DOUBLE PRECISION P,V
      SAVE /PYJETS/
	

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,XY,W2,XQ2,U
      SAVE /LEPTOU/
      REAL CUT,PARL,X,XY,W2,XQ2,U

      COMMON /ARHIDE/ PHAR(400),MHAR(400)
      SAVE /ARHIDE/
      REAL PHAR
#include "rgluco.inc"
#include "rgdisdif.inc"
#include "rgparas.inc"
#include "rgpara1.inc"
#include "rgrapgki.inc"
#include "rgpart.inc"
#include "rgherup.inc"


C...Step counter
      MSTA(4)=MSTA(4)+1
c      write(6,*) ' in arexec ..................... ',MSTA(1),IPRO
c      call pylist(1)

C...Reset some stuff
      MSTA(13)=0
      MHAR(121)=0
      MHAR(129)=0
      MHAR(135)=0
      MHAR(136)=0
      MHAR(139)=0
      NHQ=0

C...Error if ARINIT has not been called
      IF (MSTA(2).EQ.0) CALL ARERRM('AREXEC',12,0)

C...Unmark unused special positions in dipole record
      DO 100 I=MAXPAR-4,MAXPAR-2
        QQ(I)=.FALSE.
 100  CONTINUE

C...If ariadne mode just pass event through to ARPARS
      IF (MSTA(1).EQ.0) THEN
        CALL ARPARS(1,N)

C...If JETSET mode should work by just passing event on to ARPARS
      ELSEIF (MSTA(1).EQ.1) THEN
        CALL ARPARS(1,N)

C...If JETSET Matrix element mode call ARMEPS
      ELSEIF (MSTA(1).EQ.4) THEN
        CALL AREEME

C...If PYTHIA mode tag extended partons etc.
      ELSEIF (MSTA(1).EQ.2) THEN

        CALL ARPYTH

C...If LEPTO mode tag extended partons
      ELSEIF (MSTA(1).EQ.3) THEN
        IF (MSTA(32).LT.0) THEN
          CALL ARILDC
          GOTO 900
        ELSEIF (MSTA(32).GT.0) THEN
          CALL ARLEPT
          GOTO 900
        ENDIF

C...Boost to hadronic cm to avoid precision problems
        CALL ARBOLE(THEL,PHI1,PHI2,DBXL,DBYL,DBZL)

        IF (LST(24).EQ.1) THEN

          IF (MSTA(30).LT.2) THEN
            K(5,4)=0
          ELSE
            K(5,4)=3
            PARA(13)=SQRT(XQ2)
          ENDIF
          IF (MSTA(30).EQ.0) THEN
            K(6,4)=1
          ELSE
            K(6,4)=2
            PARA(12)=PARA(11)/(1.0-X)
          ENDIF
          CALL ARPARS(5,6)
        ELSEIF (LST(24).EQ.3) THEN

          IF (MSTA(30).LT.2) THEN
            K(5,4)=0
          ELSE
            K(5,4)=3
            PARA(13)=SQRT(XQ2)
          ENDIF
          IF (MSTA(30).EQ.0) THEN
            K(6,4)=1
          ELSE
            K(6,4)=2
            PARA(12)=PARA(11)/(1.0-X)
          ENDIF
          CALL ARPARS(5,6)
          IF (MSTA(30).LT.2) THEN
            K(7,4)=0
          ELSE
            K(7,4)=3
            PARA(13)=SQRT(XQ2)
          ENDIF
          IF (MSTA(30).EQ.0) THEN
            K(8,4)=1
          ELSE
            K(8,4)=2
            PARA(12)=PARA(11)/(1.0-X)
          ENDIF
          CALL ARPARS(7,8)
        ENDIF
        CALL PYROBO(1,N,0.0D0,PHI2,0.0D0,0.0D0,0.0D0)
        CALL PYROBO(1,N,THEL,PHI1,DBXL,DBYL,DBZL)

      ELSEIF(MSTA(1).EQ.5) THEN
         INPOM = 2
         XRA = sngl(Q2/dble(yy)/sss)
         IF(IPRO.EQ.12) THEN
c here QPM
            N1 = 0
            N2 = 0
            DO 30 I=3,N
               IF(IDIR.EQ.0.AND.IABS(K(I,2)).EQ.100) INPOM = I
               IF(K(I,1).EQ.2.AND.N1.EQ.0) THEN
                  N1 = I
                  DO 20 J=N1,N
                     IF(K(J,1).EQ.2.AND.K(J+1,1).EQ.1.AND.N2.EQ.0)
     +               THEN
                        N2 = J+1
                        GOTO 30
                     ENDIF
   20             CONTINUE
               ENDIF
   30       CONTINUE
C...Boost to hadronic cm to avoid precision problems
            DEL=P(NIA1,4)+P(INPOM,4)
            DBXL=(P(NIA1,1)+P(INPOM,1))/DEL
            DBYL=(P(NIA1,2)+P(INPOM,2))/DEL
            DBZL=(P(NIA1,3)+P(INPOM,3))/DEL
            CALL PYROBO(0,N,0.D0,0.D0,-DBXL,-DBYL,-DBZL)
            IF(N1.NE.0.AND.N2.NE.0) THEN
               IF(MSTA(30).LT.2) THEN
                  K(N1,4)=0
               ELSE
                  K(N1,4)=3
                  PARA(13)=SQRT(SNGL(Q2))
               ENDIF
               IF(MSTA(30).EQ.0) THEN
                  K(N2,4)=1
               ELSE
                  K(N2,4)=2
                  PARA(12)=PARA(11)/(1.0-XRA)
                  IF(IDIR.EQ.0) PARA(12)=PARA(11)/(1.0-XRA/XFGKI)
               ENDIF
c         write(6,*) ' arexec ipro=12,n1,n2 ',n1,n2
               CALL ARPARS(N1,N2)
            ENDIF
            CALL PYROBO(0,N,0.D0,0.D0,DBXL,DBYL,DBZL)
         ELSEIF(IPRO.EQ.13.OR.IPRO.EQ.14.OR.
     &	   IPRO.EQ.15.OR.IPRO.EQ.18) THEN
c            write(6,*) ' arexec 1st ',IPRO,n1,n2,n3,n4 
c here  BGF
            N1 = 0
            N2 = 0
            DO 50 I=3,N
               IF(IDIR.EQ.0.AND.IABS(K(I,2)).EQ.100) INPOM = I
               IF(K(I,1).EQ.2.AND.N1.EQ.0) THEN
                  N1 = I
                  DO 40 J=N1,N
                     IF(K(J,1).EQ.2.AND.K(J+1,1).EQ.1.AND.N2.EQ.0)
     +               THEN
                        N2 = J + 1
                        GOTO 60
                     ENDIF
   40             CONTINUE
               ENDIF
   50       CONTINUE
   60       CONTINUE
C...Boost to hadronic cm to avoid precision problems
            DEL=P(NIA1,4)+P(INPOM,4)
            DBXL=(P(NIA1,1)+P(INPOM,1))/DEL
            DBYL=(P(NIA1,2)+P(INPOM,2))/DEL
            DBZL=(P(NIA1,3)+P(INPOM,3))/DEL
            CALL PYROBO(0,N,0.D0,0.D0,-DBXL,-DBYL,-DBZL)
            IF(MSTA(30).EQ.0) THEN
               K(NIA2+1,4)=1
            ENDIF
            IF(N2.NE.0) THEN
               NN2 = N2
               N3 = 0
               N4 = 0
               DO 80 I=NN2,N
                  IF(K(I,1).EQ.2.AND.N3.EQ.0) THEN
                     N3 = I
                     DO 70 J=N3,N
                        IF(K(J,1).EQ.2.AND.K(J+1,1).EQ.1.AND.N4.EQ.0)
     +                  THEN
                           N4 = J+1
                           GOTO 90
                        ENDIF
   70                CONTINUE
                  ENDIF
   80          CONTINUE
            ENDIF
   90       CONTINUE
            IF(N1.NE.0.AND.N2.NE.0) THEN
               IF(MSTA(30).LT.2) THEN
                  K(N1,4)=0
                  IF(IDIR.EQ.0.AND.(N1+1).LT.N2) K(N2,4)=0
               ELSE
                  K(N1,4)=3
                  IF(IDIR.EQ.0.AND.(N1+1).LT.N2) K(N2,4)=3
                  PARA(13)=SQRT(SNGL(Q2))
               ENDIF
               IF(MSTA(30).EQ.0) THEN
                  K(N1+1,4)=1
               ELSE
                  K(N1+1,4)=2
                  PARA(12)=PARA(11)/(1.0-XRA)
                  IF(IDIR.EQ.0) PARA(12)=PARA(11)/(1.0-XRA/XFGKI)
               ENDIF
               CALL ARPARS(N1,N2)
            ENDIF

c         write(6,*) ' arexec 2nd ipro=13,14  ,n3,n4 ',n3,n4
            IF(N3.NE.0.AND.N4.NE.0) THEN
               IF(MSTA(30).LT.2) THEN
                  K(N3,4)=0
               ELSE
                  K(N3,4)=3
                  PARA(13)=SQRT(SNGL(Q2))
               ENDIF
               IF(MSTA(30).EQ.0) THEN
                  K(N4,4)=1
               ELSE
                  K(N4,4)=2
                  PARA(12)=PARA(11)/(1.0-XRA)
                  IF(IDIR.EQ.0) PARA(12)=PARA(11)/(1.0-XRA/XFGKI)
               ENDIF
               CALL ARPARS(N3,N4)
            ENDIF

            CALL PYROBO(0,N,0.D0,0.D0,DBXL,DBYL,DBZL)
         ENDIF


      ENDIF

C...Perform fragmentation if requested
 900  IF (MHAR(145).NE.0) CALL ARQQ2O
      IF (MSTA(5).EQ.1.and.ILHA.eq.0) CALL PYEXEC

      RETURN

C**** END OF AREXEC ****************************************************
      END
