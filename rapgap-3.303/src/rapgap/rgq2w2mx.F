C-------------------------------------------------------------------------
      SUBROUTINE RGQ2W2MX(X)
*-- Author :   H. Kowalski
C     Routine generates variables (Q2,W2,Mx2) and the
C     corresponding Jacobian factor AJAC in Common /CQ2W2MX/
C     The variables (Q2,W2,Mx2) are generated between limits
c     given in COMMON /GDLIMIT/
C     The routine computes also from Q2,W2,Mx2 the  variables
c     YBj, XBj, Xpom and Beta and puts them in COMMON /GDVARB1/
c
c     it generates also the variable Tdf (which plays a role of t)
c     Tdf generation is preliminary because it is to primitive
C-----------------------------------------------------------------------
      Implicit NONE
      Double Precision X(20),XG(20)
      Integer NDIMC
      COMMON /DIMEN/ NDIMC

      Double Precision Q2mins, Q2maxs, W2min, W2max, Mxmin, Mxmax
      COMMON /GDLIMIT/ Q2mins, Q2maxs, W2min, W2max, Mxmin, Mxmax

      Double Precision  Stot, Q2s, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2s, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision  R02, LAM, SIGM0, X0, ALPHAS
      COMMON /SATURM/   R02, LAM, SIGM0, X0, ALPHAS

      DOUBLE PRECISION Q2JAC, Yjac, MxJAC

      Double Precision  RANUNI
      External          RANUNI
      Double Precision  TWOPI
      Double Precision  W, Mx, Rsta, Rend, Tmin, Tend, mprot

      Double Precision RNMX,RNT
      Double Precision YX,FGAM,FWEI,W12,W02
      DOUBLE PRECISION ME,MP
      Integer NRN,I,J
#include "rgluco.inc"
#include "rgpara1.inc"
#include "rglujets.inc"
#include "rgpart.inc"
#include "rghs45.inc"
#include "rgdiffr.inc"
#include "rgpdiss.inc"
      Double Precision draprn
      Double precision pi

      Double Precision PHI
      COMMON/DIFFA/ PHI
	Double Precision DEMINP
      LOGICAL DEBUG
      PI = 4.D0*DATAN(1.D0)
      DEBUG = .FALSE.
C...  GIVE BEAM  FOUR VECTORS
      DO 10 I=1,2
         DO 10 J=1,5
            K(I,J)=KBEAM(I,J)
   10 P(I,J) = PBEAM(I,J)
      N=2
      ME =0.511e-3
      MP =0.938
      TWOPI = 6.2831853
c     write(6,*) ' RGQ2W2MX ',IHERAC
      NRN = 0

      IF(IHERAC.EQ.0) THEN


C
C- generate the Q2 variable
C-
         NRN = 0
         NRN = NRN + 1
         Ybj = ymin*(Ymax/Ymin)**X(NRN)
C
C- generate the Ybj variable
C-
         Ybj = ymin*(Ymax/Ymin)**X(NRN)
         YJAC = Ybj * Dlog(Ymax/Ymin)
C-
         NRN = NRN + 1
         W02=(1.+MP)**2
         W12=W02-MP*MP
         Q2MIN=ME*ME*Ybj*Ybj/(1.D0-Ybj)
         Q2MAX=Ybj*SSS - W12
         IF(QMI.GT.Q2MIN) Q2MIN = QMI
         IF(QMA.LT.Q2MAX) Q2MAX = QMA
         Q2mins = Q2MIN
         Q2maxs = Q2MAX
         Q2=Q2MIN*((Q2MAX/Q2MIN)**X(NRN))
         Q2s = Q2
         Q2JAC = Q2*DLOG(Q2MAX/Q2MIN)
         Xbj = Q2/Ybj/Stot
         PHI = 2.D0*PI*draprn()
	   nrn = nrn + 2
	   IF (NFRAG.EQ.10 .OR. NFRAG.EQ.20) NRN = NRN + 1
         NDIMC = NRN
         CALL draprnV(XG,5)
         DO I=1,5
            X(NDIMC+I)=XG(I)
         ENDDO
      ELSEIF(IHERAC.EQ.1) THEN
         NRN = 2
         Ybj = YHS
         XBj = XHS
         Q2 = Q2HS
         Q2s = Q2
         YJAC = 1.
         Q2JAC = 1.

      ELSE
         write(6,*) ' RGQ2W2Mx: wrong process IHERAC = ',IHERAC
      ENDIF

      YX = Ybj
      W2 = -Q2 + ybj*Stot
      W=sqrt(W2)

      XMAX=1.d0-XF
      if(XBJ.GE.XMAX) THEN
         AJAC = 0.
         RETURN
      ENDIF
      Mxmin = 0.3

      Mxmax = max(Mxmin,dsqrt(Q2*(XMAX-XBJ)/XBJ))
c     write(6,*) ' rgQ2W2MX: Mxmin Mxmax ',Mxmin,Mxmax,xbj,stot
c     write(6,*) ' rgQ2W2MX: ymin,ymax ',ymin,ymax,ybj,q2

      CALL PARTI(KE,YX,FGAM,FWEI,1,0)

c     write(6,*) ' RGQ2W2MX '
c     call dulist(1)
c     write(6,*) ' rgq2w2mx ',FGAM,FWEI,yx,q2
c     write(6,*) ' rgq2w2mx min ',ymin,ymax,q2min,q2max

      IF(IHERAC.EQ.0) THEN
         YJAC = YJAC * FGAM
      ENDIF
C-
C- generate the Mx variable
C-

      NRN = NRN + 1
      RNMX = DLOG(Mxmin) + X(NRN)*(DLOG(Mxmax) - DLOG(Mxmin))
      Mx = DEXP(RNMX)
      MxJAC = Mx * DLOG(Mxmax/Mxmin)

      Mx2 = Mx**2

      AJAC = Q2jac*Yjac*Mxjac

c  Generate the (positive) t variable in the range t=Tmin to t= Tend
c      Tend = 5.
      Tend = T2MAX
      If(Tend.Gt.Q2) Tend = Q2 - 0.1
      Tend = -Tend
      Tmin = 0.
      Rsta = DEXP(Tend)
      Rend = DEXP(Tmin)
      NRN = NRN + 1
      RNT = Rsta + X(NRN)*(Rend-Rsta)
      Tdf = -DLOG(RNT)/6.
c      IF (NFRAG.EQ.10 .OR. NFRAG.EQ.20) Tdf = -DLOG(RNT)/2.
cccc      NDIMC = NRN
c       write(6,*) 'rgq2w2mx MX,t ',Mx2,tdf
c       write(6,*) 'rgq2w2mx x(3),x(4) ',(x(i),i=1,ndimc)
c       Tdf = 0.
C-
C-  compute the variables xb, xg, eta
C-
      Xbj = Q2/STOT/YBj
      beta = Q2/(Mx2 + Q2 + Tdf)
      Xpom = (Mx2 + Q2 + Tdf)/(W2 + Q2)
C-
      R02 = (XPOM/X0)**LAM
c      write(*,1001) Q2,Mx,W,R02,Tdf
10000 Format('Q2 ,Mx, W, R02', 5f10.2)
c  Tmin condition
      mprot = 0.9383
      tmin = mprot**2*Xpom**2/(1.-Xpom)
	
      IF (NFRAG.EQ.10 .OR. NFRAG.EQ.20) THEN
cy MYMIN taken from gdibeg.f in diffvm:
cy DEMINP = SNGL (DMN + DMPI0 - DMP) + 0.10 
cy DEMINP = min. energy released in decay of N
cy assume dissociated proton is pion + p to get mymin
cy or could as well set mymin>mp
cccc        
         DEMINP=PYMASS(2112)+PYMASS(111)-PYMASS(2212)+ 0.10 
         MYMIN=Mprot+DEMINP
         MYMAX=20.
         MY2MIN=MYMIN*MYMIN
         MY2MAX=MYMAX*MYMAX
         NRN = NRN + 1
         MY2=MY2MIN*((MY2MAX/MY2MIN)**X(NRN))
         TMIN=-mprot**2*xpom+(MY2*xpom)/(1.D0-xpom)
c         write(6,*) 'rgq2w2mx:MY2MIN,MY2,MY2MAX :', MY2MIN, MY2, MY2MAX
c	   write(6,*) X(NRN)
c	   write(6,*) ' peps ',peps,prn2,tdf
         Ajac=Ajac * MY2**(peps) /DLOG(MY2MAX/MY2MIN)
         Ajac=Ajac * dexp(-prn2*dabs(Tdf))/dexp(-6*dabs(Tdf))
           
      ELSE
c elastic case 
         tmin = mprot**2*Xpom**2/(1.-Xpom)
      ENDIF

      If(Tdf.lt.tmin) Then
         AJAC = 0.
c         write(*,*) '***  t < Tmin  ***', Tdf, tmin
      Endif
C-

   20 RETURN
      END
