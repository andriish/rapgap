      SUBROUTINE DJGSPH(L0,IMODE)
C***********************************************************************
C
C  DJGSPH
C
C  H.Spiesberger (05.07.99)
C  Interface to SOPHIA: Simulation of events with low mass and low Q2
C  for SOPHIA 1.5
C  update 08.08.2005 for djangoh, version 1.6
C
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSCHNN/ ICHNN
      COMMON /DJPASS/ NTOT,NPASS,NFAILL
      COMMON /SPPASS/ NSOPH,NSPOUT,NFAILP,NSPACC
      COMMON /DJFAIL/ NFAILI(10)
#include "hepevt.inc"
c      PARAMETER (NMXHEP=2000)
c      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
c     &                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),
c     &                PHEP(5,NMXHEP),VHKK(4,NMXHEP)
      COMMON /S_PLIST/ P_S(2000,5), LLIST(2000), NP, Ideb
#include "rglujets.inc"
c      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
c      REAL*4 P,V
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30)
     &,XSCH,YSCH,W2SCH,Q2SCH,USCH
      REAL*4 CUT,PARL,XSCH,YSCH,W2SCH,Q2SCH,USCH
      SAVE
      DIMENSION P_gam(4),P_nuc(4)
      LOGICAL FIRST
      DATA FIRST/.TRUE./

C...Header
      IF (FIRST) THEN
        FIRST=.FALSE.
        WRITE(6,9)
 9      FORMAT(
     1'**************************************************',
     2'*****************************',
     3//,10X,'                      Calling SOPHIA '
     4//,10X,' Event generator for photohadronic processes in'
     4,' astrophysics '
     6//,10X,'                       VERSION 1.5, '/
     7//,10X,'          by A.M"ucke, Ralph Engel, J.P.Rachen, '
     8//,10X,'            R.J.Protheroe, and Todor Stanev'
     9//,10X,'         Comput.Phys.Commun. 124 (2000) 290-314',//
     9' **************************************************',
     1'****************************',//)
      CALL INITIAL(13)
      ENDIF

C...Event accepted by HERACLES
      NTOT=NTOT+1
C...reset counters for failed trials
      DO 1 INF=1,10
 1    NFAILI(INF)=0
      NSPACC=0

C...Determine momentum of virtual incoming photon
      IF (ICHNN.LE.2) THEN
C...non-radiative channels
        P_gam(1)=    -PHEP(1,1)
        P_gam(2)=    -PHEP(2,1)
        P_gam(3)=PELE-PHEP(3,1)
        P_gam(4)=EELE-PHEP(4,1)
        Q2=-T
      ELSE
C...radiative channels
        P_gam(1)=    -PHEP(1,1)-PHEP(1,3)
        P_gam(2)=    -PHEP(2,1)-PHEP(2,3)
        P_gam(3)=PELE-PHEP(3,1)-PHEP(3,3)
        P_gam(4)=EELE-PHEP(4,1)-PHEP(4,3)
        Q2=-TS
      ENDIF

C...Determine momentum of incoming nucleon
      P_nuc(1)=0D0
      P_nuc(2)=0D0
      P_nuc(3)=-PPRO
      P_nuc(4)=EPRO

      N=0
C...Generate hadronic final state
      LST(21)=99
      write(6,*) ' before sophia '
      CALL SOPHIA(L0,P_gam,P_nuc,Imode)
      call print_event(1)
      write(6,*) ' after sohia '
      pause
C...Event passed fragmentation
      IF (NP.NE.0) THEN
        LST(21)=0
        LST(36)=2
        NSOPH=NSOPH+1
        NSPACC=1
      ENDIF
      IF (NP.EQ.0) NFAILP=NFAILP+1

C...Transfer event to JETSET common block
      IF (LLEPT.EQ.-1) THEN
        LEPIN=11
      ELSEIF (LLEPT.EQ.1) THEN
        LEPIN=-11
      ENDIF
      IP1=0
C...Set initial state, first lepton:
      IP1=IP1+1
      K(IP1,1)=21
      K(IP1,2)=LEPIN
      K(IP1,3)=0
      K(IP1,4)=0
      K(IP1,5)=0
      P(IP1,1)=0.
      P(IP1,2)=0.
      P(IP1,3)=PELE
      P(IP1,4)=EELE
      P(IP1,5)=PYMASS(LEPIN)
C...then nucleon = proton
      IP1=IP1+1
      K(IP1,1)=21
      K(IP1,2)=2212
      K(IP1,3)=0
      K(IP1,4)=0
      K(IP1,5)=0
      P(IP1,1)=0.
      P(IP1,2)=0.
      P(IP1,3)=-PPRO
      P(IP1,5)=PYMASS(K(IP1,2))
      P(IP1,4)=EPRO

C...Virtual boson (non-radiative case)
      IP1=IP1+1
      K(IP1,1)=21
      K(IP1,2)=23
      K(IP1,3)=1
      K(IP1,4)=0
      K(IP1,5)=0
      IF(ICHNN.EQ.1.OR.ICHNN.EQ.2) THEN
        DO 10 JG=1,4
          P(IP1,JG)=P(1,JG)-PHEP(JG,1)
 10     CONTINUE
        P(IP1,5)=-SQRT(Q2)
      ELSE
C...Virtual boson (radiative case)
        DO 11 JG=1,4
          P(IP1,JG)=P(1,JG)-PHEP(JG,1)-P(JG,3)
 11     CONTINUE
        P(IP1,5)=-SQRT(Q2)
      ENDIF

C...scattered lepton
      IP1=IP1+1
      DO 12 J=1,5
        P(IP1,J)=PHEP(J,1)
 12   CONTINUE
      K(IP1,1)=ISTHEP(1)
      K(IP1,2)=IDHEP(1)
      K(IP1,3)=1
      K(IP1,4)=JDAHEP(1,1)
      K(IP1,5)=JDAHEP(2,1)
C...radiative photon 
      IF (ICHNN.GT.2) THEN
        IP1=IP1+1
        DO 13 J=1,5
          P(IP1,J)=PHEP(J,3)
 13     CONTINUE
        K(IP1,1)=ISTHEP(3)
        K(IP1,2)=IDHEP(3)
        K(IP1,3)=1
        IF(ICHNN.EQ.6.OR.ICHNN.EQ.8.OR.ICHNN.EQ.12) THEN
          K(IP1,3)=1
        ELSE
          K(IP1,3)=4
        ENDIF
        K(IP1,4)=JDAHEP(1,3)
        K(IP1,5)=JDAHEP(2,3)
      ENDIF

C...Hadronic final state from SOPHIA
      DO 20 IP=1,NP
        IF (ABS(LLIST(IP)).LT.50) THEN
          IP1=IP1+1
          P(IP1,1)=P_S(IP,1)
          P(IP1,2)=P_S(IP,2)
          P(IP1,3)=P_S(IP,3)
          P(IP1,4)=P_S(IP,4)
          P(IP1,5)=P_S(IP,5)
          K(IP1,1)=1
          K(IP1,2)=ICON_SIB_PDG(LLIST(IP))
          K(IP1,3)=0
          K(IP1,4)=0
          K(IP1,5)=0
        ENDIF
 20   CONTINUE
      N=IP1
      RETURN

      END
#include "datdec.f"
#include "param_ini.f"
