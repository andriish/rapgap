CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C---INTERFACE FOR CALLS OF PARTON DISTRIBUTIONS FROM HERACLES
C---TAKEN FROM PYSTFU (4.1.91) HS
C    changed for the use with pomeron structure function RASTFU
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION FUNX(XINT)
	Implicit None
      DOUBLE PRECISION T2MIN,T2,XR,MP,XP2T,XX,WTDIST,XMAX
      REAL XPQ(-6:6)
	Double Precision draprn,PYMASS
      DOUBLE PRECISION WTTEST,T2MN,T2MAX1,WTPO
      DOUBLE PRECISION XMAXV,XMINV,VMIN,VMAX,SMALL
#include "rgluco.inc"
#include "rgdiffr.inc"
#include "rgrapgki.inc"
#include "rgpdiss.inc"
	Integer KPART,IVM
      REAL SCALE
      COMMON/PINT/ SCALE,KPART
      DOUBLE PRECISION POM,PIM
      COMMON/WEIGHT1/ POM,PIM
      COMMON/VMESON/IVM
	Real FUNX,WGMAX,XINT
	Integer I,KF
      EXTERNAL PYMASS
      DATA MP/0.938/
      DATA SMALL/1.d-3/
      XMAX=1.d0-XF
      FUNX = 0.
C ... XP2 = E_glue /E_proton
C ... XX  = E_glue/E_pomeron
C ... XR  = E_pomeron/E_proton
ctest      XR = DBLE(XINT)
      XR = DBLE(XPR)*(XMAX/DBLE(XPR))**DBLE(XINT)
      IF(IVM.GE.1) THEN
         IF(IVM.GE.1.AND.IVM.LT.443) THEN
            VMAX = 1.020 + 2.d0*PYMASS(211)
            VMIN = 0.780 - 2.d0*PYMASS(211)
         ELSEIF(IVM.EQ.443) THEN
            VMAX = 2.D0*PYMASS(421)
            VMIN = 2.D0*PYMASS(4)
         ELSEIF(IVM.EQ.553) THEN
            VMAX = 2.D0*PYMASS(521)
            VMIN = 2.D0*PYMASS(5)
         ENDIF
         VMIN = VMIN+SMALL
         VMAX = VMAX-SMALL
c                  write(6,*) 'pomstr :',XMAX,XMIN
         XMAXV = DBLE(XPR)*(1 + VMAX**2/DBLE(SCALE))
         XMINV = DBLE(XPR)*(1 + VMIN**2/DBLE(SCALE))
c                  write(6,*) 'pomstr:',XMAXV,XMINV
         XR = XMINV*(XMAXV/XMINV)**DBLE(XINT)
         IF(XMAXV.GT.XMAX) RETURN
      ENDIF
      XP2T = DBLE(XPR)
      XX = XP2T/XR
      FUNX = 0.
      IF(XX.GE.0.99999) RETURN
      T2 = -99999.D0
c elastic case 
      T2MIN=MP*MP*XR*XR/(1.d0-XR)
      T2MN = -(DBLE(SCALE) *(1.D0-XR/DBLE(XPR))
     +         +4.D0*PYMASS(211)**2)

      IF(T2MN.LT.T2MAX) THEN
         T2MAX1 = T2MN
      ELSE
         T2MAX1 = T2MAX
      ENDIF

      T2 = T2MIN*((T2MAX1/T2MIN)**draprn())
c      write(6,*) ' FUNX: XR,XPR,XX',XR,XPR,XX
      IF(T2MIN.GE.T2MAX1) RETURN
      IF(XP2T.GE.XMAX) RETURN
      IF(XP2T.GE.XR) RETURN
c      write(6,*) ' pomstr ',NG,NPOM
      IF(NG.EQ.20.AND.NPOM.EQ.20) THEN
         KF=211
      ELSEIF(NG.EQ.21.AND.NPOM.EQ.21) THEN
         KF=111
      ELSE
         KF=100
      ENDIF
c      if(KF.EQ.211) THEN
c        write(6,*) ' FUNX: xx,XPR,XR ',XX,XPR,XR
c        write(6,*) ' FUNX: T2min = ',T2MIN,'T2MAX ',T2MAX
c        FUNX = 0
c        RETURN
c        endif
c      write(6,*) ' FUNX: XR,XP2T,T2mIN,T2',XR,XP2T,T2MIN,T2
c this is called to get maximum of differential pomeron distribution
      CALL RAT2DI(KF,XR,-T2,WTTEST)
      CALL RAT2IN(KF,XR,XP2T,T2MIN,-T2,WTDIST)
      IF(XX.GE.1.) THEN
         write(6,*) ' funx error ',xx,xr,xp2t
      endif
      CALL RASTFU(KF,SNGL(XX),SCALE,XPQ)
c      write(6,*) xpq(kpart),xpq(-kpart),kpart,SNGL(XX)
      FUNX = 0.
c       write(6,*) ' FUNX: WTDIST',WTDIST

      FUNX = XPQ(KPART)*SNGL(WTDIST)
      IF(IVM.EQ.0) THEN
         WTPO = T2*DLOG(T2MAX/T2MIN)*XR*DLOG(XMAX/XP2T)*WTTEST
         FUNX = FUNX * SNGL(XR) * ALOG(SNGL(XMAX)/XPR)
      ELSE
         FUNX = FUNX * SNGL(XR * DLOG(XMAXV/XMINV))
         WTPO = T2*DLOG(T2MAX/T2MIN)*XR*DLOG(XMAXV/XMINV)*WTTEST
      ENDIF
      WGMAX=0.
      DO 10  I=-6,6
         IF(I.EQ.0) GOTO 10
         WGMAX = MAX(WGMAX,XPQ(I))
   10 CONTINUE
      WTPO=WTPO*DBLE(WGMAX)
      IF(WTPO.GT.POM) THEN
         POM = WTPO
ccc         write(6,*) ' WTPO,xr,t2 ',WTPO,XR,T2
      ENDIF
      RETURN

      END
      FUNCTION FXT2(XINT,T2INT)
	Implicit None
      DOUBLE PRECISION T2MIN,T2,XR,MP,XP2T,XX,WTDIST,XMAX,WTPI
      REAL XPQ(-6:6),WGMAX
      DOUBLE PRECISION XMAXV,XMINV,VMIN,VMAX,SMALL
#include "rgluco.inc"
#include "rgdiffr.inc"
#include "rgrapgki.inc"
#include "rgpdiss.inc"
	Double Precision DEMINP
	Real SCALE
	Integer KPART,IVM
      COMMON/PINT/ SCALE,KPART
      DOUBLE PRECISION POM,PIM
      COMMON/WEIGHT1/ POM,PIM
      COMMON/VMESON/IVM
	Integer I,KF
	Double Precision T2MAX1,T2MN
	Real FXT2,XINT,T2INT
	Double Precision draprn,PYMASS
      EXTERNAL PYMASS,draprn
      DATA MP/0.938/
      DATA SMALL/1.d-3/
      XMAX=1.d0-XF
C ... XP2 = E_glue /E_proton
C ... XX  = E_glue/E_pomeron
C ... XR  = E_pomeron/E_proton
c      XR = DBLE(XINT)
      FXT2 = 1.
c      write(6,*) NG,NPOM
      IF(NG.EQ.20.AND.NPOM.EQ.20) THEN
         KF=211
      ELSEIF(NG.EQ.21.AND.NPOM.EQ.21) THEN
         KF=111
      ELSE
         KF=100
      ENDIF
      XR = DBLE(XPR)*(XMAX/DBLE(XPR))**DBLE(XINT)
      FXT2 = SNGL(XR * DLOG(XMAX/DBLE(XPR)))
      IF(IVM.GE.1) THEN
         IF(IVM.GE.1.AND.IVM.LT.443) THEN
            VMAX = 1.020 + 2.d0*PYMASS(211)
            VMIN = 0.780 - 2.d0*PYMASS(211)
         ELSEIF(IVM.EQ.443) THEN
            VMAX = 2.D0*PYMASS(421)
            VMIN = 2.D0*PYMASS(4)
         ELSEIF(IVM.EQ.553) THEN
            VMAX = 2.D0*PYMASS(521)
            VMIN = 2.D0*PYMASS(5)
         ENDIF
         VMIN = VMIN+SMALL
         VMAX = VMAX-SMALL
c                  write(6,*) 'pomstr:',XMAX,XMIN
         XMAXV = DBLE(XPR)*(1 + VMAX**2/DBLE(SCALE))
         XMINV = DBLE(XPR)*(1 + VMIN**2/DBLE(SCALE))
c                  write(6,*) 'pomstr:',XMAX,XMIN
         XR = XMINV*(XMAXV/XMINV)**DBLE(XINT)
         FXT2 = SNGL(XR * DLOG(XMAXV/XMINV))
         IF(XMAXV.GT.XMAX) RETURN
      ENDIF
      T2MN = -(DBLE(SCALE) *(1.D0-XR/DBLE(XPR))
     +         +4.D0*PYMASS(211)**2)
      IF(T2MN.LT.T2MAX) THEN
         T2MAX1 = T2MN
      ELSE
         T2MAX1 = T2MAX
      ENDIF
c generate MY2 for proton dissociation
      IF (NFRAG.EQ.10 .OR. NFRAG.EQ.20) THEN
cy MYMIN taken from gdibeg.f in diffvm:
cy DEMINP = SNGL (DMN + DMPI0 - DMP) + 0.10 
cy DEMINP = min. energy released in decay of N
cy assume dissociated proton is pion + p to get mymin
cy or could as well set mymin>mp
cccc        
         DEMINP=PYMASS(2112)+PYMASS(111)-PYMASS(2212)+ 0.10 
         MYMIN=MP+DEMINP
         MYMAX=20.
         MY2MIN=MYMIN*MYMIN
         MY2MAX=MYMAX*MYMAX
         MY2=MY2MIN*((MY2MAX/MY2MIN)**draprn())
         T2MIN=-MP*MP*XR+(MY2*XR)/(1.D0-XR)
c      T2MIN=MP*MP*XR*XR/(1.d0-XR)
c         T2MN = -(Q2  - XR*XG1*SSS  + 4.D0*PYMASS(211)**2)
c            write(6,*) 'partdf:MY2MIN,MY2,MY2MAX :', MY2MIN, MY2, MY2MAX
      ELSE
c elastic case 
      T2MIN=MP*MP*XR*XR/(1.d0-XR)
      ENDIF

      T2 = T2MIN*((T2MAX1/T2MIN)**DBLE(T2INT))
c      write(6,*) ' FXT2 : T2,T2MAX,T2MIN,T2INT ',T2,T2MAX,T2MIN,T2INT
c      T2 = DBLE(T2INT)
      XP2T = DBLE(XPR)
      XX = XP2T/XR

      IF(T2MIN.GE.T2MAX1) THEN
c          write(6,*) ' T2MIN = ',T2MIN,' T2MAX ',T2MAX
         FXT2 = 0.0
         RETURN
      ENDIF
      IF(XP2T.GE.XMAX) THEN
         FXT2 = 0.0
         RETURN
      ENDIF

      if(XR.Le.0.0) THEN
         write(6,*) ' FXT2: xx,XPR,XR ',XX,XPR,XR
         write(6,*) ' FXT2: T2min = ',T2MIN,'T2MAX ',T2MAX
c        FUNX = 0
c        RETURN
      endif
      T2GKI = SNGL(-T2)
      XFGKI = SNGL(XR)
      CALL RAT2DI(KF,XR,-T2,WTDIST)
      CALL RASTFU(KF,SNGL(XX),SCALE,XPQ)

c      write(6,*) ' FXT2: WTDIST =',wtdist
c      write(6,*) 'FXT2 ',xpq(kpart),xpq(-kpart),kpart,xx,scale
      IF (NFRAG.EQ.10 .OR. NFRAG.EQ.20) THEN
c	   write(6,*) ' peps ',peps,prn2
         WTDIST=WTDIST * MY2**(peps) /DLOG(MY2MAX/MY2MIN)
cpt      write(6,*)'weight after my distribution',wtdist
c do the t reweighting here, because wghtt is only filled after a rastfu call	   
         if(wghtt.gt.0) then
c            write(6,*) wghtt
           WTDIST = WTDIST * dexp(-prn2*dabs(t2))/wghtt
	   endif  
      ENDIF
      FXT2= FXT2 * SNGL(WTDIST)
      FXT2 = FXT2 * SNGL(T2 * DLOG(T2MAX1/T2MIN))
      IF(IVM.EQ.0) THEN
         WTPI = WTDIST*T2*DLOG(T2MAX1/T2MIN)*XR*DLOG(XMAX/XP2T)
      ELSE
         WTPI = WTDIST*T2*DLOG(T2MAX1/T2MIN)*XR*DLOG(XMAXV/XMINV)
      ENDIF
ctest      WTPI = WTDIST*T2*XR
      IF((NG.GE.30.AND.NG.LT.100).OR.(NG.LT.0.AND.NPOM.LT.0)) THEN
         WGMAX=0.
         DO 10  I=-6,6
            IF(I.EQ.0) GOTO 10
            WGMAX = MAX(WGMAX,XPQ(I))
   10    CONTINUE
         WTPI = WTPI*DBLE(WGMAX)
         IF(WTPI.GT.POM) POM=WTPI
      ELSE
         IF(WTPI.GT.PIM) PIM=WTPI
      ENDIF
      FXT2 = FXT2 * XPQ(KPART)

c      write(6,*) ' FXT2: FXT2 = ',FXT2
      RETURN
      END
      FUNCTION FLOW(XR)
      FLOW = 0.
      RETURN
      END
      FUNCTION FUP(XR)
      FUP = 1.
      RETURN
      END
      SUBROUTINE RAT2IN(KF,XR,XP2,T2MIN,T2,WTDIST)
	Implicit None
c give t distribution in proton
C KF = 100 for pomeron
C KF = 211 for pi
C XR = fractional energy of pomeron
C XP2 = fractional energy of parton to hard interaction
C T2 = t = (p - p')**2 < 0
C WTDIST = weight of distribution
#include "rgdiffr.inc"

#include "rgpara.inc"
      INTEGER KF
      DOUBLE PRECISION T2,XR,XP2,T2MIN,WTDIST,XE,XMAX,T2CH,beta02
      LOGICAL FIRST

      DATA FIRST/.TRUE./
      DATA beta02/58.74D0/
      XMAX=1.d0-XF

      WTDIST = 0.D0
c      write(6,*) 'RAT2IN: XR,XP2,T2MIN,T2',XR,XP2,T2MIN,T2
      IF(NPOM.EQ.0) THEN
         IF(FIRST) THEN
            write(6,*) ' streng pomeron is used'
            FIRST = .FALSE.
         ENDIF
c pomeron distribution from streng (hera proc. 1987)
         XE=2.D0*ALPHP*DLOG(1.d0/XR)+RN2
         T2CH = T2MAX
         IF(XE*T2MIN.GT.170.D0) THEN
c                      WRITE(6,*) 'XE..',XE,DEXP(-XE)
            GOTO 10
         ENDIF
         IF(XE*T2MAX.GT.170.D0) T2CH = 100.D0/XE
c         write(6,*) ' XE ',xe,' XR ',XR
c         write(6,*) 'T2MAX ',T2MAX,' T2MIN ',T2MIN
         WTDIST = DEXP(-XE*T2MIN)-DEXP(-XE*T2CH)
c         write(6,*) WTDIST
         WTDIST = WTDIST/XE
         WTDIST = WTDIST*beta02/XR**(1.D0 + 2.D0 * EPSP) /16.d0/PI
         IF(DABS(WTDIST).LT.1.D-20) GOTO 10
      ELSEIF(NPOM.EQ.1) THEN
         IF(FIRST) THEN
            write(6,*) ' Ingelman pomeron is used'
            FIRST = .FALSE.
         ENDIF
c test with Ingelman f_pom distribution
         T2CH = T2MAX
         IF(8.D0*T2MAX.GT.170.D0) T2CH = 170.D0/8.D0
         WTDIST = 6.38D0/8.D0*(DEXP(-8.D0*T2MIN) - DEXP(-8.D0*T2CH))
         WTDIST = WTDIST+ 0.424D0/8.D0*(DEXP(-3*T2MIN) - DEXP(-3*T2CH))
         WTDIST = WTDIST/XR/2.3D0
c end Ingelman
      ELSE
         write(6,*) ' RAT2IN: pomeron distribution ',NPOM,
     +       ' not implemented'
         write(6,*) ' RAT2IN: program stops '
         stop
      ENDIF
      RETURN
   10 WTDIST = 0.D0
      RETURN
      END
