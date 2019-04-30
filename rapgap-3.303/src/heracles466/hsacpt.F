C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSACPT(ICONT)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSKPXY/ XX,Y
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSCHNN/ ICHNN
      COMMON /HSONLY/ IHSONL
#include "hepevt.inc"
c      PARAMETER (NMXHEP=2000)
c      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
c     &                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),
c     &                PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      COMMON /SOPHCT/ WSOPHIA
      COMMON /SPPASS/ NSOPH,NSPOUT,NFAILP,NSPACC
      LOGICAL IELAST
      LOGICAL IDJSPH
      DIMENSION IQFLAV(-6:6)
      DIMENSION IQFLCC(-6:6)
      DATA IQFLAV /-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6/
      DATA IQFLCC /-5,-6,-3,-4,-2,-1,0,1,2,4,3,6,5/
      DATA NCEVE /0/

      IELAST=.FALSE.
      IDJSPH=.FALSE.
      ICHNN=ICONT
      IF (ICONT.EQ.3.OR.ICONT.EQ.15.OR.ICONT.EQ.16.OR.ICONT.EQ.17)
     &  IELAST=.TRUE.
      X=XX

C---------------------------------------      COUNT ACCEPTED EVENTS
      NCEVE=NCEVE + 1
      NEVHEP=NCEVE

C----------------------------------------     HADRONIC MASS
      GSP=SP-MEI2-MPRO2
      IF (ICONT.GT.5) THEN
        OMEGA=(PH*DKQ+PQH*DKP)/(PH*EQH+PQH*EH)
        DCTHGA=(DKP/OMEGA-ME2/2D0/EH)/PH
        CTHGA=1D0-DCTHGA
        DKPRO=OMEGA*(EPRO+CTHGA*PPRO)
      ENDIF
      IF (ICONT.LE.2) THEN
        W2=Y*(1D0-X)*GSP+MPRO2
      ELSEIF (IELAST) THEN
        W2=MPRO2
      ELSE
        W2=Y*(1D0-X)*GS+MPRO2-2D0*(DKP-DKPS+DKPRO)
      ENDIF
      W=DSQRT(W2)
      IF (W.LE.WSOPHIA) IDJSPH=.TRUE.
      IFAIL=0
      IF ((.NOT.IDJSPH).AND.IHSONL.EQ.0) CALL HSWCUT(W,1,12,IFAIL)
      IF (IFAIL.NE.0) RETURN

C----------------------------------------     SAMPLE QUARK TYPE
      IF (IPDFOP.GE.1.AND.(.NOT.IELAST).AND.(.NOT.IDJSPH)) THEN
  100   CALL HSFLAV(SNGL(W2),IQF,IQFR)
C...Check if energy in jet system is enough for fragmentation.
        IFAIL=0
        IF (IHSONL.EQ.0) CALL HSWCUT(W,IQF,IQFR,IFAIL)
        IF (IFAIL.NE.0) GOTO 100
        ELSE
        IQF=0
        IQFR=0
      ENDIF

C...test print
      IF(NEVHEP.LE.10.AND.IPRINT.GE.3) THEN
        WRITE(LUNTES,'(///A/5X,2I5,1PE15.4)')
     *        ' HSACPT : NEVHEP, IQF, CQP(12)',
     *                   NEVHEP, IQF
        WRITE(LUNTES,'(2(5X,6(1PE15.4)/))') CQP
      ENDIF

C------------------------------         RECONSTRUCT KINEMATICS
      IF(ICONT.LE.3) THEN
C------------------------------         NON-RADIATIVE CHANNELS
        Q2=X*Y*GSP
C -----------------------------         DEFINE  COMMON /HSIKP/
        S=X*SP
        T=-Q2
        U=-X*SP+Q2
        SS=S
        TS=T
        US=U
        DKP=0D0
        DKPS=0D0
        DKQ=0D0
        DKQS=0D0
C
        AMFEL=(SP-MPRO2-MEI2)/2D0-Q2/2D0/X
        BMFEL=Q2/2D0+MEI2
        EFEL=(BMFEL*PPRO+AMFEL*PELE)/(PELE*EPRO+PPRO*EELE)
        PFEL=SQRT((EFEL-ME)*(EFEL+ME))
        CTFEL=EELE/PELE*EFEL/PFEL - BMFEL/PELE/PFEL
        STFEL2=1D0-CTFEL*CTFEL
        IF (STFEL2.GT.0D0) THEN
          STFEL=DSQRT(STFEL2)
          ELSE
          STFEL=0D0
        ENDIF
        CALL SFECFE(SPHIEL,CPHIEL)
        EFQU=EELE+X*EPRO-EFEL
        PFQU=EFQU
        CTFQU=(PELE-X*PPRO-PFEL*CTFEL)/PFQU
        STFQU=SQRT((1D0+CTFQU)*(1D0-CTFQU))
C
C------------------------------        FILL STANDARD COMMON
        NHEP=4
C------------------------------        FINAL ELECTRON : IHEP=1
        IHEP=1
        ISTHEP(IHEP)=1
        IF (ICONT.EQ.1.OR.ICONT.EQ.3) THEN
          IDHEP(IHEP)=-LLEPT*11
          ELSE
          IDHEP(IHEP)=-LLEPT*12
        ENDIF
        PHEP(1,IHEP)=CPHIEL*STFEL*PFEL
        PHEP(2,IHEP)=SPHIEL*STFEL*PFEL
        PHEP(3,IHEP)=CTFEL*PFEL
        PHEP(4,IHEP)=EFEL
        PHEP(5,IHEP)=ME
        JMOHEP(1,IHEP)=0
        JMOHEP(2,IHEP)=0
        JDAHEP(1,IHEP)=0
        JDAHEP(2,IHEP)=0
C-------------------------------       HADRONIC FINAL STATE
C-------------------------------       FOR (QUASI)ELASTIC EP:
C                                      IT IS A PROTON
        IF (IELAST) THEN
          IHEP=2
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=2212
          PHEP(1,IHEP)=-CPHIEL*STFQU*PFQU
          PHEP(2,IHEP)=-SPHIEL*STFQU*PFQU
          PHEP(3,IHEP)=CTFQU*PFQU
          PHEP(4,IHEP)=EFQU
          PHEP(5,IHEP)=0D0
          JMOHEP(1,IHEP)=0
          JMOHEP(2,IHEP)=0
          JDAHEP(1,IHEP)=0
          JDAHEP(2,IHEP)=0
C-------------------------------       IF PARTON MODEL IS USED:
C                                      SCATTERED QUARK
        ELSEIF (IPDFOP.GE.1) THEN
          IHEP=2
          ISTHEP(IHEP)=1
          IF (ICONT.EQ.1) THEN
            IDHEP(IHEP)=IQFLAV(IQF)
            ELSE
            IDHEP(IHEP)=IQFLCC(IQF)
          ENDIF
          PHEP(1,IHEP)=-CPHIEL*STFQU*PFQU
          PHEP(2,IHEP)=-SPHIEL*STFQU*PFQU
          PHEP(3,IHEP)=CTFQU*PFQU
          PHEP(4,IHEP)=EFQU
          PHEP(5,IHEP)=0D0
          JMOHEP(1,IHEP)=0
          JMOHEP(2,IHEP)=0
          JDAHEP(1,IHEP)=0
          JDAHEP(2,IHEP)=0
C-------------------------------       SPECTATOR QUARK SYSTEM
          IHEP=4
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=90
          PHEP(1,IHEP)=          -(PHEP(1,1)+PHEP(1,2))
          PHEP(2,IHEP)=          -(PHEP(2,1)+PHEP(2,2))
          PHEP(3,IHEP)=-PPRO+PELE-(PHEP(3,1)+PHEP(3,2))
          PHEP(4,IHEP)= EPRO+EELE-(PHEP(4,1)+PHEP(4,2))
          PHEP(5,IHEP)=0D0
        ELSE
C-------------------------------    HADRONIC FINAL STATE : IHEP=2
C                                   FOR THE STRUCTURE FUNCTIONS OPTION
          IHEP=2
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=93
          PHEP(1,IHEP)=-CPHIEL*STFQU*PFQU
          PHEP(2,IHEP)=-SPHIEL*STFQU*PFQU
          PHEP(3,IHEP)=PELE-PPRO-PHEP(3,1)
          PHEP(4,IHEP)=EELE+PPRO-PHEP(4,1)
          PHEP(5,IHEP)=W
          JMOHEP(1,IHEP)=0
          JMOHEP(2,IHEP)=0
          JDAHEP(1,IHEP)=0
          JDAHEP(2,IHEP)=0
C-------------------------------       PROTON REMNANT NOT SEPARATED
          IHEP=4
          ISTHEP(4)=0
          DO 1 IAX=1,5
    1     PHEP(IAX,IHEP)=0D0
        ENDIF
C-------------------------------       NO PHOTON : IHEP=3
        IHEP=3
        ISTHEP(IHEP)=0
        DO 2 IAX=1,5
    2   PHEP(IAX,IHEP)=0D0
C
C------------------------------------------------------------------------------
C------------------------------   ***  HARD PHOTON CONTRIBUTIONS
      ELSEIF(ICONT.GT.5) THEN
        Q2=X*Y*GSP
        NHEP=4
C-------------------------------       FINAL ELECTRON : IHEP=1
        IHEP=1
        ISTHEP(IHEP)=1
        IF (ICONT.LE.10.OR.ICONT.GE.15) THEN
          IDHEP(IHEP)=-LLEPT*11
          ELSE
          IDHEP(IHEP)=-LLEPT*12
        ENDIF
        PHEP(4,IHEP)=ESH
        PHEP(3,IHEP)=COSEH*PSH
        PHEP(1,IHEP)=SINEH*PSH
        PHEP(2,IHEP)=0D0
        PHEP(5,IHEP)=ME
C-------------------------------       FINAL STATE PHOTON : IHEP=3
        IHEP=3
        ISTHEP(IHEP)=1
        IDHEP(IHEP)=22
C---ALREADY DONE:
C       OMEGA=(PH*DKQ + PQH*DKP) / (PH*EQH + PQH*EH)
C       DCTHGA=(DKP/OMEGA-ME2/2D0/EH)/PH
C       CTHGA=1D0-DCTHGA
        STHGA=DSQRT(DCTHGA*(2D0-DCTHGA))
        PPP2=PSH*SINEH
        STCPHG=(ESH - DKPS/OMEGA - PSH*CTHGA*COSEH) / PPP2
        PHEP(4,IHEP)=OMEGA
        PHEP(3,IHEP)=OMEGA*CTHGA
        PHEP(1,IHEP)=OMEGA*STCPHG
        AKY2=PHEP(4,IHEP)**2-PHEP(1,IHEP)**2-PHEP(3,IHEP)**2
        IF (AKY2.LT.0D0) THEN
          PHEP(2,IHEP)=0D0
          ELSE
          PHEP(2,IHEP)=DSQRT(AKY2)*DSIGN(1D0,0.5D0-HSRNDM(-1))
        ENDIF
        PHEP(5,IHEP)=0D0
C-------------------------------       HADRONIC FINAL STATE
        IF (IELAST) THEN
C-------------------------------
C                                      FOR (QUASI)ELASTIC EP: PROTON
          IHEP=2
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=2212
          ETOT=EELE+EQH
          PZTOT=PELE-PQH
          PHEP(4,IHEP)=ETOT-PHEP(4,1)-PHEP(4,3)
          PHEP(1,IHEP)=-PHEP(1,1)-PHEP(1,3)
          PHEP(2,IHEP)=-PHEP(2,1)-PHEP(2,3)
          PHEP(3,IHEP)=PZTOT-PHEP(3,1)-PHEP(3,3)
          PHEP(5,IHEP)=0D0
        ELSEIF (IPDFOP.GE.1) THEN
C-------------------------------
C                                      FOR PARTON MODEL:
C                                      SCATTERED QUARK : IHEP=2
C                                      USE E-P CONSERVATION WITHOUT
C                                      QUARK PT
          IHEP=2
          ISTHEP(IHEP)=1
          IF (ICONT.LE.10) THEN
            IDHEP(IHEP)=IQFLAV(IQF)
            ELSE
            IDHEP(IHEP)=IQFLCC(IQF)
          ENDIF
          ETOT=EELE + EQH
          PZTOT=PELE - PQH
          PHEP(4,IHEP)=ETOT - PHEP(4,1) - PHEP(4,3)
          PHEP(1,IHEP)= - PHEP(1,1) - PHEP(1,3)
          PHEP(2,IHEP)= - PHEP(2,1) - PHEP(2,3)
          PHEP(3,IHEP)=PZTOT - PHEP(3,1) - PHEP(3,3)
          PHEP(5,IHEP)=0D0
C-------------------------------       SPECTATOR QUARK SYSTEM
          IHEP=4
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=90
          PHEP(1,IHEP)=-(PHEP(1,1)+PHEP(1,2)+PHEP(1,3))
          PHEP(2,IHEP)=-(PHEP(2,1)+PHEP(2,2)+PHEP(2,3))
          PHEP(3,IHEP)=-PPRO+PELE
     &                 -(PHEP(3,1)+PHEP(3,2)+PHEP(3,3))
          PHEP(4,IHEP)= EPRO+EELE
     &                 -(PHEP(4,1)+PHEP(4,2)+PHEP(4,3))
          PHEP(5,IHEP)=0D0
        ELSE
C-------------------------------    HADRONIC FINAL STATE : IHEP=2
C                                   FOR STRUCTURE FUNCTIONS OPTION
          IHEP=2
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=93
          ETOT=EELE+EPRO
          PZTOT=PELE-PPRO
          PHEP(4,IHEP)=ETOT-PHEP(4,1)-PHEP(4,3)
          PHEP(1,IHEP)=-PHEP(1,1)-PHEP(1,3)
          PHEP(2,IHEP)=-PHEP(2,1)-PHEP(2,3)
          PHEP(3,IHEP)=PZTOT-PHEP(3,1)-PHEP(3,3)
          PHEP(5,IHEP)=W
C-------------------------------       PROTON REMNANT NOT SEPARATED
          IHEP=4
          ISTHEP(4)=0
          DO 3 IAX=1,5
    3     PHEP(IAX,IHEP)=0D0
        ENDIF
C-------------------------------       ROTATE WHOLE EVENT
        PHIEVE=HSRNDM(-1)*2D0*PI
        CROT=DCOS(PHIEVE)
        SROT=DSIN(PHIEVE)
        DO 4 INHEP=1,4
        PNHEPX= CROT*PHEP(1,INHEP)+SROT*PHEP(2,INHEP)
        PNHEPY=-SROT*PHEP(1,INHEP)+CROT*PHEP(2,INHEP)
        PHEP(1,INHEP)=PNHEPX
    4   PHEP(2,INHEP)=PNHEPY
C
      ELSE
        STOP
      ENDIF
C
      IF (IHSONL.EQ.0) THEN
        CALL DJGVAR(ICHNN,X,Y,Q2)
        NSPACC=0
        IF (IDJSPH) THEN
          CALL DJGSPH(13,IMODE)
        ELSE
          CALL DJGEVT
        ENDIF
      ENDIF
C-------------------------
C---Define beam particles for /HEPEVT/
chs...ISTHEP according to H1 standard
      NHEP = NHEP + 2
C---incoming electron:
      IHEP=NHEP-1
      ISTHEP(IHEP)=201
      IDHEP(IHEP)=-LLEPT*11
      PHEP(5,IHEP)=0D0
      PHEP(4,IHEP)=EELE
      PHEP(1,IHEP)=0D0
      PHEP(2,IHEP)=0D0
      PHEP(3,IHEP)=PELE
      VHEP(4,IHEP)=0D0
      VHEP(1,IHEP)=0D0
      VHEP(2,IHEP)=0D0
      VHEP(3,IHEP)=0D0
      JMOHEP(1,IHEP)=0
      JMOHEP(2,IHEP)=0
      JDAHEP(1,IHEP)=0
      JDAHEP(2,IHEP)=0
C incoming proton:
      IHEP = NHEP
      ISTHEP(IHEP)=201
      IDHEP(IHEP)=2212
      PHEP(5,IHEP)=0D0
      PHEP(4,IHEP)=EPRO
      PHEP(1,IHEP)=0D0
      PHEP(2,IHEP)=0D0
      PHEP(3,IHEP)=-PPRO
      VHEP(4,IHEP)=0D0
      VHEP(1,IHEP)=0D0
      VHEP(2,IHEP)=0D0
      VHEP(3,IHEP)=0D0
      JMOHEP(1,IHEP)=0
      JMOHEP(2,IHEP)=0
      JDAHEP(1,IHEP)=0
      JDAHEP(2,IHEP)=0
C
      CALL HSUSER(2,X,Y,Q2)
C
      RETURN
      END
