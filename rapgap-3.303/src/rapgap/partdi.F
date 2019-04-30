      SUBROUTINE PARTDI(X,WPART)
      IMPLICIT NONE
#include "rgrapgki.inc"
#include "rglujets.inc"
#include "rgluco.inc"
#include "rgpart.inc"
#include "rgpara1.inc"
#include "rgpara.inc"
#include "rgdiffr.inc"
#include "rgstru.inc"
#include "rghs45.inc"
#include "rgnlo.inc"
      Double Precision PHI
      COMMON/DIFFA/ PHI
      Double Precision draprn
      DOUBLE PRECISION ME,MP
      Double Precision X,XY,WPART,WEIGHT
      COMMON/XVAR/ XY(10)
      Double Precision  DTHETA,DPHI
      DIMENSION X(20)
      REAL XPQ(-6:6),XPQF2(-6:6)
      REAL F2CC(-6:6),XF3CC(-6:6)
      Double Precision WMAX
      Integer NDIMC,IMIX,IDEBUG,NDIM,NPOIN,LST,IRES,NTHRE,NVMRE,IGENFL
      COMMON /DIMEN/ NDIMC
      COMMON /OALPHAS/ WMAX,IMIX
      COMMON/ INTERN/IDEBUG
      COMMON/DIVO/ NDIM,NPOIN
      COMMON/EPPARA/LST(30),IRES(2)
      COMMON/HERTHET/ NTHRE,NVMRE
      COMMON/GENWEI/IGENFL
      DOUBLE PRECISION QG2
      COMMON/SEMIH/ QG2
      REAL SNGL
      EXTERNAL draprn
      Double Precision QF,QFT,ALPH_EM,GF,DOT,WTGLU
      Integer NFLAVP,NACC
      Double Precision  WT,W02,W12,YX,SMIN,XP2Q,XG1,XP1MIN
      Double Precision QG2MAX,QG2MIN,QG20,XZOHRN,WQG2
      Double Precision  XP2,SHAT1,PIO,FGAM,FWEI,FLUX
      Double Precision THETE,ECM,PT2,WTG15,XP1,PR,t,at
      Double Precision SPHP,CPHP,STHP,CTHP,BOCHCK
      Double Precision WMAT,XPINT,XP,XX,QQ
	Double Precision F2BGF
	External F2BGF
      Integer I,J,IN,KI,NP,NPP,KPA2,KPAT,KPFL,KPAO,KPAO2
      Integer NDIMEN,NRN,IST,NPFIN,NB2,KPF,KFLCH,KFLSP
      Integer KPZ
      

      DATA QF/0.0D0/,W12/0.0D0/
      NDIMEN = NDIM
c      IDEBUG=1
      SCAL1 = -99999.
      XPD1 = -99999.
      SCAL2 = -99999.
      XPD2 = -99999.
      DO 10  IN=1,20
         K(IN,1) = 0
         K(IN,2) = 0
   10 CONTINUE

C...  GIVE BEAM  FOUR VECTORS
      DO 20 I=1,2
         DO 20 J=1,5
            K(I,J)=KBEAM(I,J)
   20 P(I,J) = PBEAM(I,J)
      ME =P(1,5)
      MP =P(2,5)
      NP=2
      NPP=3*NP
      N=2

      QF = 0.D0
      NFLAVP = 0

      IF(IPRO.EQ.10.OR.IPRO.EQ.12.OR.IPRO.EQ.13.
     +OR.IPRO.EQ.15.OR.IPRO.EQ.18) THEN
         AM(1) = 0.0D0
         AM(2) = 0.0D0
         IF(IHFLA.GE.4.AND.IPRO.EQ.18) THEN
            KPA=IHFLA
            AM(1) = PYMASS(KPA)
            AM(2) = PYMASS(KPA)
         ENDIF
      ELSEIF(IPRO.EQ.16) THEN
         AM(1) = 0.0D0
         AM(2) = 0.0D0
      ELSEIF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
         IF(IWEI.EQ.0) THEN
            KPA = 4
            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
         ENDIF
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA)
      ELSEIF(IPRO.EQ.102) THEN
	   KPA = 25
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.500) THEN
         KPZ = 23
         KPA = KPZ
         AM(1) = PYMASS(KPA)
         AM(2) = 0.
      ELSEIF(IPRO.EQ.501) THEN
         KPZ = 23
         KPA = KPZ
         KPA2= 1
         AM(1) = PYMASS(KPA)
         AM(2) = 0.d0
      ELSEIF(IPRO.EQ.503) THEN
         NP=3
         NPP=3*NP
	   KPA = 1
         KPZ = 23
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA)
         AM(3) = PYMASS(KPZ)
      ELSEIF(IPRO.EQ.504) THEN
         NP=3
         NPP=3*NP
	   KPA = IHFLA
         KPZ = 23
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA)
         AM(3) = PYMASS(KPZ)
      ELSEIF(IPRO.EQ.510) THEN
      ELSEIF(IPRO.EQ.513) THEN
         NP=3
         NPP=3*NP
c        select charge of W
         KPZ = 24
	   KPA = -2
	   KPA2 = 1
         if(draprn().ge.0.5) then 
            KPZ=-24         
	      KPA = 2
            KPA2 = -1
         Endif   
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA2)
         AM(3) = PYMASS(KPZ)
      ELSEIF(IPRO.EQ.514) THEN
         NP=3
         NPP=3*NP
c        select charge of W
         KPZ = 24
         if(IHFLA.eq.5) then
            write(6,*) ' p_semih: cabbibo supressed combination selected'
            write(6,*) ' p_semih: heaviest flavor = ',IHFLA
            write(6,*) ' p_semih: not implemented - Program stopped '
            stop
            endif
	   KPA = -IHFLA
         KPA2 = IHFLA-1
         if(draprn().ge.0.5) then 
            KPZ=-24         
	      KPA = IHFLA
            KPA2 = -(IHFLA-1)
         Endif   
         AM(1) = PYMASS(KPA)
         AM(2) = PYMASS(KPA2)
         AM(3) = PYMASS(KPZ)
      ELSE
         WRITE(6,*) 'wrong subprocess selected: IPRO = ',IPRO
         WRITE(6,*) '**** PROGRAM STOP ****'
         STOP
      ENDIF
      IF(IPRO.EQ.14.AND.IMIX.EQ.1) THEN
      ELSE
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
         W02=(AM(1)+AM(2)+MP)**2
         IF(IPRO.EQ.12) W02=(AM(1)+MP)**2
         IF(AM(1).LT.1.0D0) W02=(1.D0 + MP)**2
         IF(AM(1).LT.1.0D0.and.pt2cut(ipro).ne.0.)
     +     W02=(DSQRT(PT2CUT(IPRO)) + MP)**2
         W12=W02-MP*MP
         YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
         YMAX=YMAX/(2.D0*(SSS+ME*ME))
         YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
         YMIN=YMIN/(2.D0*(SSS+ME*ME))
         IF(YMI.GT.YMIN) YMIN=YMI
         IF(YMA.LT.YMAX) YMAX=YMA
      ENDIF
c         write(6,*) w02
c         WRITE(6,10500) YMIN,YMAX
c10500 FORMAT(' limits on y ',/,
c     +' YMIN = ',E10.5,' YMAX = ',E10.5)

C ... select particle code for light flavour production according
C ... to charges
      IF(IPRO.EQ.10.OR.IPRO.EQ.13.OR.IPRO.EQ.18) THEN
         IF(IHFLA.GE.4.AND.IPRO.EQ.18) THEN
            IF(IWEI.EQ.0) THEN
               KPA=4
               IF(IHFLA.GE.4) KPA=IHFLA
            ELSEIF(IWEI.EQ.1) THEN
               IF(draprn().GT.0.5) KPA = -KPA
            ENDIF
         ELSEIF(INTER.LT.2) THEN
            QF=DFLOAT(PYCHGE(1))**2 + DFLOAT(PYCHGE(2))**2 + DFLOAT(
     +      PYCHGE(3))**2
            QF = 2.D0*QF
            KPA = -4
            QFT = - draprn()*QF
   30       KPA=KPA+1
            QFT = QFT + DBLE(PYCHGE(KPA))**2
            IF(QFT.LT.0.0D0) GOTO 30
            IF(KPA.GT.3) write(6,*) 'fatal light quark = charm!!!!!! ',
     +      KPA
            KPA2 = - KPA
         ELSEIF(INTER.EQ.2) THEN
            KPAT = 2
            IF(draprn().GT.0.5) THEN
               KPA = -ISIGN(1,K(1,2))*IABS(KPAT)
               KPA2 = ISIGN(1,K(1,2))*(IABS(KPAT)-1)
            ELSE
               KPA = ISIGN(1,K(1,2))*(IABS(KPAT)-1)
               KPA2 = -ISIGN(1,K(1,2))*IABS(KPAT)
            ENDIF
         ENDIF
      ELSEIF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
         IF(IWEI.EQ.0) THEN
            KPA=4
            IF(IHFLA.GE.4) KPA=IHFLA
         ELSEIF(IWEI.EQ.1) THEN
            IF(draprn().GT.0.5) KPA = -KPA
         ENDIF
         KPFL = IABS(KPA)
         KPA2=-KPA
         IF(INTER.EQ.2) THEN
            KPA = -ISIGN(1,K(1,2))*KPFL
            KPA2 = ISIGN(1,K(1,2))*(KPFL-1)
            IF(draprn().LT.0.5) THEN
               KPA = ISIGN(1,K(1,2))*(KPFL-1)
               KPA2 = -ISIGN(1,K(1,2))*KPFL
            ENDIF
         ENDIF
         AM(1)=PYMASS(KPA)
         AM(2)=PYMASS(KPA2)
      ELSEIF(IPRO.EQ.12.OR.IPRO.EQ.15.OR.IPRO.EQ.16) THEN
         NFLAVP = NFLAV
         KPA = 1
      ELSEIF(IPRO.eq.501) THEN
      ELSEIF(IPRO.ge.503.and.IPRO.le.505) THEN
         IF(IWEI.EQ.0) THEN

         ELSEIF(IWEI.EQ.1) THEN
            IF(draprn().GT.0.5) KPA = -KPA
         ENDIF
         KPA2=-KPA
         AM(2)=PYMASS(KPA2)
         AM(1)=PYMASS(KPA)
      ELSEIF(IPRO.ge.513.and.IPRO.le.515) THEN
         AM(2)=PYMASS(KPA2)
         AM(1)=PYMASS(KPA)
      ENDIF
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
      ELSEIF(IGENFL.EQ.1) THEN
         KPA = KPAO
         KPA2 = KPAO2
         IF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
            AM(1)=PYMASS(KPA)
            AM(2)=PYMASS(KPA2)
         ENDIF
      ENDIF

C... YX IS THE PHOTON ENERGY
C... Q2 FROM PHOTON
C... XP2 IS XGLUON ( MOMENTUM FRACTION OF THE GLUON)
C... XMIN = MIN XGLUON TO PRODUCE THE INV. MASS OF GAMMA GLUON SYSTEM
C... XMAX=1.
CCC
C... GENERATE YX,Q2,XP2 ACCORDING TO 1/X SPECTRUM
C... FGAM IS THE WEIGHT OF EPA
      NRN = 0
      XMAX=0.999d0
      IF(IRES(1).EQ.0.OR.IPRO.EQ.12) THEN
c.......
         IF(KE.NE.22) THEN
            NRN = NRN + 1

            YX = YMIN*((YMAX/YMIN)**X(NRN))
            IF(YX.GT.YMAX.OR.YX.LT.YMIN) GOTO 200
            Q2MIN=ME*ME*YX*YX/(1.D0-YX)
            IF(QMI.GT.Q2MIN) Q2MIN = QMI
            Q2MAX=YX*SSS - W12
            IF(QMA.LT.Q2MAX) Q2MAX = QMA
            IF(Q2MAX.LT.Q2MIN) GOTO 210
            NRN = NRN + 1
            IF(QMI.EQ.0.D0) THEN
               Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
            ELSE
c new try weighting with 1/q**4
               Q2 = Q2MIN*Q2MAX/(X(NRN)*Q2MIN + Q2MAX*(1.D0 - X(NRN)))
            ENDIF
            XMIN = 0.0D0
            XP2Q=Q2/YX/SSS
            QG2 = 0.D0
            WQG2 = 1.D0
            IF(ISEMIH.EQ.1) THEN
C... here is new for ZOTOV SALEEV
               QG2MAX = 90.D0
               IF(AM(1).GT.1.D0) THEN
                  QG2MAX = 4.D0 * AM(1)**2 + Q2
               ELSE
                  QG2MAX = 4.D0 * PT2CUT(IPRO) + Q2
               ENDIF
               QG2MIN = 0.2D0
               QG20 = 2.D0
               NRN = NRN + 1
               XZOHRN = X(NRN)
c 1/qt**2
               QG2 = QG2MIN*((QG2MAX/QG2MIN)**XZOHRN)
               WQG2 = QG2*DLOG(QG2MAX/QG2MIN)
               IF(QG2.LE.QG20) QG2 = 0.D0
c... end new ZOTOV
            ENDIF
            IF(IPRO.NE.12) THEN
               SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
               XMIN=(SMIN+Q2+QG2)/(YX*SSS)
               IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
                  XMIN = DMAX1(XMIN,Q2/YX/SSS)
               ENDIF
            ELSE
               XMIN=Q2/YX/SSS
            ENDIF
         ELSEIF(KE.EQ.22) THEN
            Q2 = 0.D0
            YX = 1.D0
            QG2 = 0.D0
            WQG2 = 1.D0
            IF(ISEMIH.EQ.1) THEN
C... here is new for ZOTOV SALEEV
               QG2MAX = 90.D0
               IF(AM(1).GT.1.D0) THEN
                  QG2MAX = 4.D0 * AM(1)**2
               ELSE
                  QG2MAX = 4.D0 * PT2CUT(IPRO)
               ENDIF
               QG2MIN = 0.2D0
               QG20 = 2.D0
               NRN = NRN + 1
               XZOHRN = X(NRN)
c 1/qt**2
               QG2 = QG2MIN*((QG2MAX/QG2MIN)**XZOHRN)
               WQG2 = QG2*DLOG(QG2MAX/QG2MIN)
               IF(QG2.LE.QG20) QG2 = 0.D0
c... end new ZOTOV
            ENDIF
            SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
            XMIN= (SMIN+QG2)/(YX*SSS)
         ELSE
            write(6,*) ' wrong KF selected; program stop '
            YX = 0.D0
         ENDIF
         XG1=YX
      ELSEIF(IRES(1).EQ.1.AND.IABS(KE).NE.2212) THEN
C XG1 = E_gluon/E_electron
C YX  = E_photon/E_electron
C XP1 = E_gluon/E_photon
         NRN = NRN + 1
         YX = YMIN*((YMAX/YMIN)**X(NRN))
         NRN = NRN + 1
         XP1MIN =DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
     +           /YX/SSS
         XP1 = XP1MIN*(1.D0/XP1MIN)**X(NRN)
         XG1 = XP1*YX
         Q2MIN = ME*ME*YX*YX/(1.D0 - YX)
         IF(QMI.GT.Q2MIN) Q2MIN = QMI
         Q2MAX=YX*(SSS-MP**2) - W12
         IF(QMA.LT.Q2MAX) Q2MAX = QMA
         IF(Q2MAX.LT.Q2MIN) GOTO 210
         NRN = NRN + 1
         IF(QMI.EQ.0.D0) THEN
            Q2 = Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
         ELSE
c new try weighting with 1/q**4
            Q2 = Q2MIN*Q2MAX/(X(NRN)*Q2MIN + Q2MAX*(1.D0 - X(NRN)))
         ENDIF
         XMIN = 0.0D0
         XP2Q=Q2/XG1/SSS
         QG2 = 0.D0
         WQG2 = 1.D0
         IF(IPRO.NE.12) THEN
            SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
            XMIN=(SMIN+Q2+QG2)/(XG1*SSS)
            IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
               XMIN = DMAX1(XMIN,Q2/YX/SSS)
            ENDIF
         ELSE
            XMIN=Q2/YX/SSS
         ENDIF
      ELSEIF(IRES(1).EQ.1.AND.IABS(KE).EQ.2212) THEN
C XG1 = E_gluon/E_proton
C YX  = XG1
C XP1 = 1
         NRN = NRN + 1
         YX = YMIN*((YMAX/YMIN)**X(NRN))
         XP1 = 1.d0
         XG1 = XP1*YX
         XP1MIN=1.d0
         Q2 = 0.d0
         XMIN = 0.0D0
         XP2Q=Q2/XG1/SSS
         QG2 = 0.D0
         WQG2 = 1.D0
         SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
         if(IPRO.ge.501.and.IPRO.lt.503) then
           SMIN = DMAX1(4.D0*PT2CUT(IPRO),(AM(1)+AM(2))**2)         
           Elseif(IPRO.ge.503.and.ipro.le.515) then
           SMIN = DMAX1(4.D0*PT2CUT(IPRO),(AM(1)+AM(2)+AM(3))**2)
           Elseif(IPRO.eq.102.or.IPRO.eq.500) then
           SMIN = AM(1)**2
           XMIN = SMIN/(XG1*SSS)
         Endif  
         XMIN=DMIN1(SMIN/(XG1*SSS),1d0)
      ELSE
         YX = 0.D0
         XG1 = 0.D0
      ENDIF
      IF(XMIN.GE.XMAX) GOTO 220
C ... XP2 = E_parton /E_proton

      IF(IPRO.eq.102.or.IPRO.EQ.500) then
         XP2 = XMIN
      ELSEIF(IPRO.NE.12.AND.IMIX.EQ.0) THEN
         NRN = NRN + 1
         XP2= XMIN*((XMAX/XMIN)**X(NRN))
      ELSEIF(IPRO.NE.12.AND.IMIX.EQ.1.AND.IWEI.EQ.1) THEN
         XP2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
      ELSEIF(IPRO.EQ.12.OR.IPRO.EQ.100) THEN
         XP2 = Q2/YX/SSS
      ELSE
         XP2 = 0.D0
      ENDIF


		
      IF(SNGL(XP2).GE.0.9999.OR.XP2.LT.XMIN) GOTO 230
      SHAT1 = SSS*XG1*XP2
      IF(IPRO.NE.12)  SHAT1=SSS*XG1*XP2 - Q2 - QG2
      IF(SHAT1.LT.(AM(1)+AM(2))**2) GOTO240
      ENTRY PARTDIHS(X,WPART)
      IF(IHERAC.EQ.1) THEN
         DO 40     IN=1,20
            K(IN,1) = 0
            K(IN,2) = 0
   40    CONTINUE

C...  GIVE BEAM  FOUR VECTORS
         DO 50    I=1,2
            DO 50     J=1,5
               K(I,J)=KBEAM(I,J)
   50    P(I,J) = PBEAM(I,J)
         ME =P(1,5)
         MP =P(2,5)
         NP=2
         NPP=3*NP
         N=2
         NFLAVP = NFLAV
         NRN = 2
         MP =P(2,5)
         IRES(1)=0
         XMAX= 0.999d0
         IST = 0
         YX = YHS
         XP2 = XHS
         Q2 = Q2HS
         IF(IPRO.EQ.12.OR.IPRO.EQ.15) THEN
            KPA = 1
         ENDIF
         XG1 = YX
         NDIMEN = NDIM
C ... select particle code for light flavour production according
C ... to charges
         IF(IPRO.EQ.10.OR.IPRO.EQ.13) THEN
            IF(INTER.LT.2) THEN
               QF=DFLOAT(PYCHGE(1))**2 + DFLOAT(PYCHGE(2))**2 +
     +         DFLOAT( PYCHGE(3))**2
               QF = 2.D0*QF
               KPA = -4
               QFT = - draprn()*QF
   60          KPA=KPA+1
               QFT = QFT + DBLE(PYCHGE(KPA))**2
               IF(QFT.LT.0.0D0) GOTO 60
               IF(KPA.GT.3) write(6,*) 'fatal light quark = charm!!!!!!'
     +         //' ', KPA
               KPA2 = - KPA
            ELSEIF(INTER.EQ.2) THEN
               KPAT = 2
               IF(draprn().GT.0.5) THEN
                  KPA = -ISIGN(1,K(1,2))*IABS(KPAT)
                  KPA2 = ISIGN(1,K(1,2))*(IABS(KPAT)-1)

               ELSE
                  KPA = ISIGN(1,K(1,2))*(IABS(KPAT)-1)
                  KPA2 = -ISIGN(1,K(1,2))*IABS(KPAT)
               ENDIF
            ENDIF
         ENDIF
         IF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
            IF(IWEI.EQ.0) THEN
               KPA = 4
               IF(IHFLA.GE.4) KPA=IHFLA
               KPA2=-KPA
            ELSEIF(IWEI.EQ.1) THEN
               IF(draprn().GT.0.5) KPA = -KPA
               KPA2 = -KPA
            ENDIF
            KPFL = IABS(KPA)
            IF(INTER.EQ.2) THEN
               KPA = -ISIGN(1,K(1,2))*KPFL
               KPA2 = ISIGN(1,K(1,2))*(KPFL-1)
               IF(draprn().LT.0.5) THEN
                  KPA = ISIGN(1,K(1,2))*(KPFL-1)
                  KPA2 = -ISIGN(1,K(1,2))*KPFL
               ENDIF
            ENDIF
            AM(1)=PYMASS(KPA)
            AM(2)=PYMASS(KPA2)
            IF(AM(1).LT.1.D-6) AM(1)=0.0D0
            IF(AM(2).LT.1.D-6) AM(2)=0.0D0
         ENDIF

         IF(IPRO.EQ.14.AND.IMIX.EQ.1) THEN
         ELSE

c         write(6,*) 'PARTDI SSS ',SSS
C.. HERE THE LIMITS ON Y( PHOTON ENERGY) AND Q**2 ARE CALCULATED
C... ACCORDING TO PROCEEDINGS OF HERA WORKSHOP 1987
C... ALI ET AL
            W02=(1.D0 + MP)**2
            W12=W02-MP*MP
            YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
            YMAX=YMAX/(2.D0*(SSS+ME*ME))
            YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
            YMIN=YMIN/(2.D0*(SSS+ME*ME))
            IF(YMI.GT.YMIN) YMIN=YMI
            IF(YMA.LT.YMAX) YMAX=YMA
         ENDIF
         IF(IGENFL.EQ.0) THEN
            KPAO = KPA
            KPAO2 = KPA2
         ELSE
            KPA=KPAO
            KPA2 = KPAO2
         ENDIF

         Q2MIN=ME*ME*YX*YX/(1.D0-YX)
         IF(QMI.GT.Q2MIN) Q2MIN = QMI
         Q2MAX=YX*SSS - W12
         IF(QMA.LT.Q2MAX) Q2MAX = QMA
c A. Kappes         IF(Q2MAX.LT.Q2MIN) GOTO 210
         XP2Q=Q2/YX/SSS
         IF(IPRO.NE.12) THEN
            SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
            QG2 = 0.D0
            IF(IPRO.EQ.100) SMIN = 2.D0*(AM(1)**2+AM(2)**2)
            XMIN=(SMIN+Q2+QG2)/(YX*SSS)
            IF(XMIN.GE.XMAX) GOTO 220
            IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
               XMIN = DMAX1(XP2Q,XMIN)
               XP2= XMIN*((XMAX/XMIN)**X(NDIMEN+1))
               if(xp2.GE.1.d0) then
                  write(6,*) ' partdi xp2,x(ndimen+1),xmax,xmin',
     +            xp2,x(ndimen+1),xmax,xmin
               endif
            ENDIF
         ENDIF
         IF(XP2.GE.XMAX) GOTO 240
      ENDIF
      IST = 0
      YY=SNGL(YX)
      XEL=SNGL(XG1)
      XPR=SNGL(XP2)
      PHI = 99999.D0
      IF(IPRO.EQ.15.OR.IPRO.EQ.16) NFLAVP = NFLQCDC
      IF(IABS(K(1,2)).NE.2212) THEN
         IF((IPRO.EQ.13.OR.IPRO.EQ.14.OR.IPRO.EQ.15.OR.IPRO.EQ.18)
     +   .AND.IMIX.EQ.0) THEN
            NRN=NRN+1
            PHI = 2.D0*PI*X(NRN)
         ELSEIF((IPRO.EQ.13.OR.IPRO.EQ.14.OR.IPRO.EQ.15.OR.IPRO.EQ.16
     +   .OR.IPRO.EQ.18)
     +   .AND.IMIX.EQ.1.AND.IWEI.EQ.1) THEN
            PHI = 2.D0*PI*X(NDIMEN+2)
         ELSEIF(IHERAC.EQ.0) THEN
            IF(IGENFL.EQ.0) THEN
               PHI = 2.D0*PI*draprn()
               PIO = PHI
            ELSEIF(IGENFL.EQ.1) THEN
               PHI = PIO
            ENDIF
         ENDIF
      ELSEIF(IABS(K(1,2)).EQ.2212) THEN
         IF(IGENFL.EQ.0) THEN
            PHI = 2.D0*PI*draprn()
            PIO = PHI
         ELSEIF(IGENFL.EQ.1) THEN
            PHI = PIO
         ENDIF
      ENDIF
      IF(IRES(1).EQ.0) THEN
         CALL PARTI(KE,YX,FGAM,FWEI,1,IST)
      ELSEIF(IRES(1).EQ.1) THEN
         CALL PARTI(KE,YX,FGAM,WEIGHT,1,IST)
c factor xg1/yx cancelled because parton densities of photon taken
c xgx instead of g(x)
         FWEI =  DLOG(1.D0/XP1MIN) * WEIGHT
         if(iabs(ke).eq.2212) FWEI = WEIGHT
      ELSE
         WRITE(6,*) ' IRES(1) > 1 not implemented '
         FGAM = 0.D0
         FWEI = 0.D0
      ENDIF
      XP=0.D0
      IRES(2)=1
C here do some gymnastics to make cut on theta angle of electron
      IF(THEMA.LT.180.D0.OR.THEMI.GE.0.D0) THEN
C go first to ep LAB system
         CALL PYROBO(3,3,0.D0,0.D0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
         THETE = PYP(3,14)
         CALL PYROBO(3,3,0.D0,0.D0,-CM(1)/CM(4),
     +                -CM(2)/CM(4),-CM(3)/CM(4))
         IF(THETE .GT. 180.01D0) THEN
            WRITE(6,*) ' FATAL: theta_electron > 180 deg ',THETE
         ENDIF
         IF(THETE .GT. THEMA) GOTO 250
         IF(THETE .LT. THEMI) GOTO 250
      ENDIF
C end of these gymnastics

C FINAL STATE PROTON
      NPFIN=NIA1+4
      N=NPFIN
      K(NPFIN,1)=1
      K(NPFIN,2)=KP
      K(NPFIN,3)=2
      P(NPFIN,5)=PYMASS(KP)
      PR = P(2,3)*(1.D0 - XP2)
      P(NPFIN,4) = PR
      CPHP=1.D0
      SPHP=DSQRT(1.D0 - CPHP**2)
      IF(ISEMIH.EQ.1) THEN
         STHP = DSQRT(QG2)/PR
      ENDIF
      IF(DABS(STHP).GT.1.D0) goto 310
      CTHP=DSQRT(1.D0 - STHP**2)
      P(NPFIN,1)= PR*STHP*CPHP
      P(NPFIN,2)= PR*STHP*SPHP
      P(NPFIN,3)= PR*CTHP
C MOMENTA OF PARTON
      K(NIA1+1,1)=21
      K(NIA1+1,2)=KGL
      IF(IPRO.EQ.12) K(NIA1+1,2) = KPA
      K(NIA1+1,3)=2
      DO 70  KI=1,4
         P(NIA1+1,KI)=P(2,KI)-P(NPFIN,KI)
   70 CONTINUE
      IF(ISEMIH.EQ.0) THEN
         P(NIA1+1,4)=ABS(P(NIA1+1,3))
         P(NIA1+1,5)=0.0D0
      ELSE
         P(NIA1+1,4)=DSQRT(P(NIA1+1,1)**2+P(NIA1+1,2)**2+P(NIA1+1,3)**2
     +               - QG2)
         P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
      ENDIF
      NIA2 = NIA1+1

      IF(IPRO.NE.12) THEN
         NF1=NIA1+2
         NF2=NIA1+3
         K(NF1,1)=2

         IF(IPRO.EQ.15) THEN

         ELSE
            K(NF2,2)=KPA2
         ENDIF
         IF(IPRO.NE.15.AND.IPRO.NE.12) KPF=KPA
         if(IPRO.ge.503.and.IPRO.le.515) NF2=NIA1+4
         K(NF1,2)=KPF
         K(NF1,3)=NIA1
         K(NF2,1)=1
         K(NF2,3)=NIA1
         P(NF1,5)=AM(1)
         P(NF2,5)=AM(2)
         if(IPRO.ge.503.and.IPRO.le.515) then
            K(NF1+1,1)=1
            K(NF1+1,2)=KPZ
            K(NF1+1,3)=NIA2
            P(NF1+1,5)=AM(3)
c            elseif(ipro.eq.15.or.ipro.eq.16) then
c            K(NF1,2)=21
c            K(NF2,2)=2 
            elseif(ipro.eq.501) then
            K(NF1,1)=1
            K(NF1,2)=KPA
            K(NF2,2)=KPA2 
         endif   
         if(kpa.eq.0) then
            write(6,*) 'p-semih kpa ',kpa
         endif
         NB2 = NIA2
C...   VECTOR OF GAMMA GLUON CMS SYSTEM
      ELSEIF(IPRO.EQ.12) THEN
         Q2Q = Q2
         CALL rgstfu(K(2,2),SNGL(XP2),SNGL(Q2Q),XPQ)
C... WTGLU IS THE WEIGHT FOR XPARTON GENERATION
         WTGLU = 0.D0
         IF(INTER.LT.2) THEN
            xx=xp2
		qq=dsqrt(q2)
            DO 80  I=-NFLAVP,NFLAVP
		XPQF2(I) = XPQ(I)
		IF(IBS.EQ.1) THEN
		   XPQF2(I) = XPQ(I) + F2BGF(I,xx,qq)
		   ENDIF
   80       WTGLU =WTGLU+DBLE(XPQF2(I))*DFLOAT(PYCHGE(I))**2/9.D0
            IF(IGENFL.EQ.1) THEN
               KPA = KPAO
               WTGLU = DBLE(XPQF2(KPA))*DFLOAT(PYCHGE(KPA))**2/9.D0
            ELSE
               QFT = - draprn()*WTGLU
               KPA=-NFLAVP-1
   90          KPA=KPA+1
               QFT = QFT + DBLE(PYCHGE(KPA))**2/9.D0*DBLE(XPQF2(KPA))
               IF(QFT.LT.0.0D0) GOTO 90
            ENDIF
            KPF = KPA
         ELSEIF(INTER.EQ.2) THEN
            DO 100 I=-6,6
               F2CC(I) = 0.0
  100       XF3CC(I) = 0.0
            DO 110 I=1,NFLAVP-1,2
               IF(ISIGN(1,K(1,2)).EQ.1) THEN
                  F2CC(I+1) = 2.*XPQ(I+1)
                  F2CC(-I) =  2.*XPQ(-I)
                  XF3CC(I+1) = 2.*XPQ(I+1)
                  XF3CC(-I) =  -2.*XPQ(-I)
                  WTGLU = WTGLU +
     +            (1.D0 - YX + YX**2/2.D0)*DBLE(F2CC(-I)+F2CC(I+1))
     +            + (YX - YX**2/2.D0)*DBLE(XF3CC(-I)+XF3CC(I+1))
               ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                  F2CC(I) = 2.*XPQ(I)
                  F2CC(-I-1) = 2.*XPQ(-I-1)
                  XF3CC(I) = 2.*XPQ(I)
                  XF3CC(-I-1) = -2.*XPQ(-I-1)
                  WTGLU = WTGLU +
     +            (1.D0 - YX + YX**2/2.D0)*DBLE(F2CC(I)+F2CC(-I-1))
     +            - (YX - YX**2/2.D0)*DBLE(XF3CC(I)+XF3CC(-I-1))
               ENDIF
  110       CONTINUE
            IF(IGENFL.EQ.1) THEN
               KPA = KPAO
               IF(ISIGN(1,K(1,2)).EQ.1) THEN
                  WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +            DBLE(F2CC(KPA)) + (YX - YX**2/2.D0)*DBLE(XF3CC(KPA))
               ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                  WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +            DBLE(F2CC(KPA)) - (YX - YX**2/2.D0)*DBLE(XF3CC(KPA))
               ENDIF
               IF(ISIGN(1,K(1,2)).EQ.1) THEN
                  KPF = (KPA - 1)
               ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                  KPF = (KPA + 1)
               ENDIF
            ELSE
               QFT = - draprn()*WTGLU
               KPA=-NFLAVP-1
  120          KPA=KPA+1
               IF(KPA.EQ.0) GOTO 120
               IF(ISIGN(1,K(1,2)).EQ.1) THEN
                  QFT = QFT + (1.D0 - YX + YX**2/2.D0)*DBLE(F2CC(KPA))
     +            + (YX - YX**2/2.D0)*DBLE(XF3CC(KPA))
               ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                  QFT = QFT + (1.D0 - YX + YX**2/2.D0)*DBLE(F2CC(KPA))
     +            - (YX - YX**2/2.D0)*DBLE(XF3CC(KPA))
               ENDIF

               IF(QFT.LT.0.0D0) GOTO 120
               IF(ISIGN(1,K(1,2)).EQ.1) THEN
                  KPF = (KPA - 1)
               ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                  KPF = (KPA + 1)
               ENDIF
            ENDIF
         ELSE
            WRITE(6,*) ' wrong interaction INTER = ',INTER
         ENDIF
         P(NIA1+1,5) = 0.0D0
         P(NIA1+1,4)=SQRT(P(NIA1+1,3)**2 + P(NIA1+1,5)**2)
         K(NIA1+1,2)=KPA
         NF1=NIA1+2
         NF2=NIA1+3
         K(NF1,1)=2
         K(NF1,2)=KPF
         K(NF1,3)=NIA1
         K(NF2,1)= 0
         NB2 = 2
         SCAL2 = SNGL(Q2Q)
         XPD2 = XPQ(KPA)
      ELSE
         NB2 = 0
      ENDIF
      DBCMS(1)=  P(NIA1,1) + P(NB2,1)
      DBCMS(2)=  P(NIA1,2) + P(NB2,2)
      DBCMS(3)=  P(NIA1,3) + P(NB2,3)
      DBCMS(4)=  P(NIA1,4) + P(NB2,4)
      DO 130 I=1,4
         P(NF1,I)=0.0D0
         P(NF2,I)=0.0D0
         if(IPRO.ge.503.and.IPRO.le.515) P(NF1+1,I)=0.0D0
  130 CONTINUE
      SHAT=DOT(DBCMS,DBCMS)
      IF(SHAT.LE.0.0) THEN
         GOTO 270
      ENDIF
C NOW BOOST TO GAMMA GLUON
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2 +
     +(DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 260
      CALL PYROBO(0,N,0.D0,0.D0,-DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +-DBCMS(3)/DBCMS(4))
      DPHI = PYANGL(P(NIA1,1),P(NIA1,2))
      call PYROBO(0,0,0.D0,-dphi,0.d0,0.d0,0.d0)
      DTHETA = PYANGL(P(NIA1,3),P(NIA1,1))
      call PYROBO(0,0,-DTHETA,0.D0,0.d0,0.d0,0.d0)
      IF(IPRO.NE.12) THEN
C NOW  LOOK THAT WE REALLY HAVE ENOUGH ENERGY IN GAMMA GLUON CMS SYSTEM
C...  ECM = CMS ENERGY OF GAMMA GLUON SYSTEM
         ECM =DSQRT(SHAT)
         If(ipro.eq.3.or.ipro.eq.102.or.ipro.eq.500)  then 
         else
           if(IPRO.ge.503.and.IPRO.le.515) then
            IF(ECM.LE.(AM(1)+AM(2)+AM(3))) GOTO 280
           else
           IF(ECM.LE.(AM(1)+AM(2))) GOTO 280
           endif
         endif
      IF(IPRO.EQ.3.or.IPRO.EQ.102.or.ipro.eq.500) Then
         K(NF2,2)=0
         K(NF2,1)=-1
         K(NF2,3)=0
         P(NF2,5)=0
         K(NF1,1)=1
         K(NF1,3)=4
        DO I=1,4
            P(NF1,I)=0
            P(NF2,I)=0
         Enddo
         P(NF1,4)=ECM
         P(NF1,5)=AM(1)
c         if(ecm/am(1).le.0.98.or.ecm/am(1).gt.1.02) then
c            write(6,*) ' partdi ecm,Am ',ecm,am(1)
c         endif
c        WT to corrrect for it afterwards	   
	   WT = (2.D0*PI)**NPP
	ELSE
         IF(IMIX.EQ.1.AND.IWEI.EQ.1) THEN
            XY(1)=X(NDIMEN+3)
            XY(2)=X(NDIMEN+4)
            if(IPRO.ge.503.and.IPRO.le.515) then
              XY(3)=X(NDIMEN+5)
              XY(4)=X(NDIMEN+6)
              XY(5)=X(NDIMEN+7)
           endif
         ELSE
            NRN = NRN + 1
            XY(1) = X(NRN)
            NRN = NRN + 1
            XY(2) = X(NRN)
            if(IPRO.ge.503.and.IPRO.le.515) then
              NRN = NRN + 1
              XY(3) = X(NRN)
              NRN=NRN+1
              XY(4)=X(NRN)
              NRN=NRN+1
              XY(5)=X(NRN)
           endif
         ENDIF
         CALL PHASE(NP,ECM,AM,PCM,WT)
         IF(WT.LE.0.D0) GOTO 290
         DO 140 I=1,4
            P(NF1,I)=PCM(I,1)
            P(NF2,I)=PCM(I,2)
            if(IPRO.ge.503.and.IPRO.le.515) then
              P(NF1+1,I)=PCM(I,3)
            endif
  140    CONTINUE
	Endif
         PT2 = PYP(NF1,9)
         CALL CUTG(PT2,NACC)
         IF(NACC.EQ.0) GOTO 300
         IF(IQ2.EQ.1) THEN
            Q2Q = (2.D0*AM(1))**2
         ELSEIF(IQ2.EQ.2) THEN
            Q2Q = SHAT
         ELSEIF(IQ2.EQ.3) THEN
            Q2Q = (2.D0*AM(1))**2 + PT2
         ELSEIF(IQ2.EQ.4) THEN
            Q2Q = Q2
         ELSEIF(IQ2.EQ.5) THEN
            Q2Q = Q2 + PT2 + (2.D0*AM(1))**2
c this is for testing
         ELSEIF(IQ2.EQ.6) THEN
c torbjorn suggestion
            Q2Q = PT2 + 8.D0*PT2**2/SHAT
C kt**2 as from Zeppenfeld/Mirkes
            t = AM(1)**2 - 2.D0 * DOT1(NIA1,NF1) - P(NIA1,5)**2
            at=DSQRT((t+Q2)**2+4.D0*Q2*PT2)
            Q2Q = at*(at+t+Q2)/2.D0/Q2
c         write(6,*) 'partdi : t,PT2,Q2,Q2Q ',t,PT2,Q2,Q2Q
         ELSEIF(IQ2.EQ.7) THEN
c torbjorn suggestion
            Q2Q = PT2 + 8.D0*PT2**2/SHAT
         ELSE
            WRITE(6,*) ' NO VALID Q2 SCALE. STOP'
            STOP
         ENDIF
         IF(IQ2.EQ.5) THEN
            Q2Q = Q2 + PT2*SCALFA + (2.D0*AM(1))**2
         ELSE
            Q2Q = SCALFA*Q2Q
         ENDIF
c now do QCD compton.....
         IF(IPRO.EQ.15.OR.IPRO.EQ.16) THEN
            IF(IPRO.EQ.15) then 
             K(NF2,2) = 21
             else
             K(NF2,2) = 22
             endif
cdon't use that            Q2Q = Q2
            if(XP2.GT.1.d0) write(6,*) 'partdi xp2>1',xp2,xmax,xmin,ipro
            CALL rgstfu(K(2,2),SNGL(XP2),SNGL(Q2Q),XPQ)
C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
            WTGLU = 0.D0
            IF(INTER.LT.2) THEN
               DO 150 I=-NFLAVP,NFLAVP
  150          WTGLU =WTGLU+DBLE(XPQ(I))*DFLOAT(PYCHGE(I))**2/9.D0
               IF(IGENFL.EQ.1) THEN
                  KPA = KPAO
                  WTGLU = DBLE(XPQ(KPA))*DFLOAT(PYCHGE(KPA))**2/9.D0
               ELSE
                  QFT = - draprn()*WTGLU
                  KPA=-NFLAVP-1
  160             KPA=KPA+1
                  QFT = QFT + DBLE(PYCHGE(KPA))**2/9.D0*DBLE(XPQ(KPA))
                  IF(QFT.LT.0.0D0) GOTO 160
               ENDIF
               KPF = KPA
            ELSEIF(INTER.EQ.2) THEN
               DO 170 I=-6,6
                  F2CC(I) = 0.0
  170          XF3CC(I) = 0.0
               DO 180 I=1,NFLAVP-1,2
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     F2CC(I+1) = 2.*XPQ(I+1)
                     F2CC(-I) =  2.*XPQ(-I)
                     XF3CC(I+1) =  2.*XPQ(I+1)
                     XF3CC(-I) = - 2.*XPQ(-I)
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(-I)+F2CC(I+1)) + (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(-I)+XF3CC(I+1))
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     F2CC(I) =  2.*XPQ(I)
                     F2CC(-I-1) =  2.*XPQ(-I-1)
                     XF3CC(I) =  2.*XPQ(I)
                     XF3CC(-I-1) = - 2.*XPQ(-I-1)
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(I)+F2CC(-I-1)) - (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(I)+XF3CC(-I-1))
                  ENDIF
  180          CONTINUE
               IF(IGENFL.EQ.1) THEN
                  KPA = KPAO
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(KPA)) + (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(KPA))
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     WTGLU = WTGLU + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(KPA)) - (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(KPA))
                  ENDIF
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     KPF = (KPA - 1)
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     KPF = (KPA + 1)
                  ENDIF
               ELSE
                  QFT = - draprn()*WTGLU
                  KPA=-NFLAVP-1
  190             KPA=KPA+1
                  IF(KPA.EQ.0) GOTO 190
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     QFT = QFT + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(KPA)) + (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(KPA))
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     QFT = QFT + (1.D0 - YX + YX**2/2.D0)*
     +               DBLE(F2CC(KPA)) - (YX - YX**2/2.D0)*
     +               DBLE(XF3CC(KPA))
                  ENDIF

                  IF(QFT.LT.0.0D0) GOTO 190
                  IF(ISIGN(1,K(1,2)).EQ.1) THEN
                     KPF = (KPA - 1)
                  ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                     KPF = (KPA + 1)
                  ENDIF
               ENDIF
            ELSE
               WRITE(6,*) ' wrong interaction INTER = ',INTER
            ENDIF
            WTG15 = WTGLU
            K(NIA1+1,2)=KPA
            AM(1) = 0.D0
            IF(IABS(KPA).GE.4) AM(1) = PYMASS(KPA)
c end QCD compton.....
         ENDIF
C BOOST BACK TO OVERALL CMS ( now done in FXN1 )
         call PYROBO(0,0,DTHETA,0.D0,0.d0,0.d0,0.d0)
         call PYROBO(0,0,0.D0,dphi,0.d0,0.d0,0.d0)
         XPR = SNGL(XP2)
         if(xp2.GT.1.d0) write(6,*) 'partdi xp2>1',xp2,ipro
         CALL rgstfu(K(2,2),SNGL(XP2),SNGL(Q2Q),XPQ)

C... WTGLU IS THE WEIGHT FOR XGLUON GENERATION
c for resolved photon move parton densities into ELERES
         IF(IRES(1).EQ.1) THEN
            XP = 1.D0/XP2
         ELSE
            XP = DBLE(XPQ(0))/XP2
         ENDIF
         SCAL2 = SNGL(Q2Q)
         XPD2 = XPQ(0)
         IF(IPRO.EQ.15.OR.IPRO.EQ.16) THEN
            XP = WTG15/XP2
            XPD2 = XPQ(KPA)
         ENDIF
         XPINT=DLOG(XMAX/XMIN)
         WTGLU = XP2 * XP * XPINT

         IF(ISEMIH.EQ.1.AND.IPRO.NE.12.AND.IPRO.NE.15.AND.IPRO.NE.16)
     +      WTGLU = WTGLU*WQG2
	   if(ipro.eq.102.or.ipro.eq.500) wtglu=1
         WMAT=WT *WTGLU*FGAM*FWEI
         WMAT=WMAT/((2.D0*PI)**NPP)
C...  FLUX = FLUX USED FOR X-SECTION CALCULATION
         FLUX=(2.D0*PI)**4
         FLUX=FLUX/(2.D0*(SHAT+Q2))
         WPART = WMAT *FLUX
         IF(IPRO.EQ.13.OR.IPRO.EQ.14.OR.IPRO.EQ.15.OR.IPRO.EQ.16) THEN
            WPART = WTGLU * FWEI * WT
            IF(IPRO.EQ.15.OR.IPRO.EQ.16) THEN
               WPART = WPART*FLUX/((2.D0*PI)**NPP)
            ENDIF

         ENDIF
	   IF(ipro.eq.102.or.ipro.eq.500) then
	      Wpart = WTGLU*FGAM*FWEI 
c xsection for resonance production gg-> h0
c sig = \int dx1 dx2 g(x1) g(x2) sigma_hat \delta(x1x2 * sss - M2)
c     = \int dx1 dx2 g(x1) g(x2) sigma_hat 1/(x1*sss) \delta(x2 - M2/(x1*sss))
c     = \int dx1 g(x1) g(x2=M2/x1*sss) sigma_hat x2/M2 
c     = \int dx1/x1  x1g(x1) x2g(x2) sigma_hat/M2
            Wpart=Wpart/AM(1)**2
         Endif
         PT2H = SNGL(PT2)
         SHH = SNGL(SHAT)

      ELSEIF(IPRO.EQ.12) THEN
         K(NIA2,2) = KPA
         K(NF1,2) = KPA
         CALL PYSPLI(K(2,2),KPA,KFLCH,KFLSP)
         IF(SHAT.LT.(PYMASS(KPF) + PYMASS(KFLSP))**2) GOTO 240
         IF(SHAT.LT.4.d0) GOTO 240
         CALL PY2ENT(NF1,KPF,KFLSP,DSQRT(SHAT))
         N = NPFIN
         K(NF1,3) = NIA1+1
         K(NF1+1,3) = 2
         K(NPFIN,1)= 0
         CALL PYEDIT(12)
         NF2 = NF1

         ALPH_EM = ALPH
         IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
         IF(INTER.LT.2) THEN
            WMAT = 2.D0 * ALPH_EM*ALPH_EM*PI/YX/Q2**2
            WMAT = WMAT *(1.D0+(1.D0-YX)**2)
         ELSEIF(INTER.EQ.2) THEN
            GF = PI*ALPH_EM/(SIN2W*XMW2*DSQRT(2.D0))
            WMAT = GF**2/2.D0/PI /(1.D0 + Q2/XMW2)**2 /YX
capply additional factor 0.5 because only right/or left handed electrons contr.
            WMAT = WMAT * 0.5D0
         ENDIF
         WPART = WMAT * WTGLU *FWEI
C BOOST BACK TO OVERALL CMS ( now done in FXN1 )
         call PYROBO(0,0,DTHETA,0.D0,0.d0,0.d0,0.d0)
         call PYROBO(0,0,0.D0,dphi,0.d0,0.d0,0.d0)
         CALL PYROBO(0,N,0.D0,0.D0,DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +   DBCMS(3)/DBCMS(4))
         PT2H = 99999.
         SHH = 99999.
      ENDIF
      NDIMC = NRN
      XFGKI = 9999.
      T2GKI = 9999.
      IF(IGENFL.EQ.0) THEN
         KPAO = KPA
         KPAO2 = KPA2
      ELSE
         KPA = KPAO
         KPA2 = KPAO2
      ENDIF
      RETURN
  200 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: ylimit ; RETURN ',yx,ymin,ymax,x(1)
      ENDIF
      NDIMC = 9999
      RETURN
  210 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: Q2limit ; RETURN ',Q2,Q2MIN,Q2MAX
         NDIMC = 9999
      ENDIF
      RETURN
  220 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: xlimit ; RETURN xmin,xmax',xmin,xmax
      ENDIF
      NDIMC = 9999
      RETURN
  230 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: xp2imit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  240 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: slimit ; RETURN ',shat1,am(1),am(2),IPRO
         write(6,*) ' sss,xg1,xp2,q2,qg2 ',sss,xg1,xp2,q2,qg2
         write(6,*) ' partdi q2,y,sss',q2,yx,sss
      ENDIF
      NDIMC = 9999
      RETURN
  250 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: theta limit ; RETURN '
      ENDIF
      NDIMC = 9999
      IF(IHERAC.EQ.1) THEN
         NTHRE = NTHRE +1
      ENDIF
      RETURN
  260 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: bochck limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  270 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: shat ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  280 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: ECM limit ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN
  290 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: PHASE WT=0 ; RETURN '
      ENDIF
      NDIMC = 9999
      RETURN

  300 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: PTCUT limit ; RETURN ',PT2,KPA
      ENDIF
      NDIMC = 9999
      RETURN
  310 CONTINUE
      WPART = 0.D0
      IF(IDEBUG.EQ.1) THEN
         write(6,*) ' partdi: QG2 limit ; RETURN ',DSQRT(QG2),PR,STHP
      ENDIF
      NDIMC = 9999
      RETURN


      END
