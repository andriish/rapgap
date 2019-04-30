      SUBROUTINE ELEQQGG(WT1)
C
C    e p  ----> Q Q_BAR GLUON GLUON p
C
C
C    full martix element
      implicit none
#include "rgpart.inc"
#include "rglujets.inc"
#include "rgpara1.inc"
#include "rgpara.inc"
#include "rgluco.inc"
#include "rgrapgki.inc"
#include "rgpqcdpom.inc"

C fuer massiv: naechste 2 Zeilen
      REAL mass
      COMMON /MASSE/mass
      DOUBLE PRECISION P0E(4),P2E(4),P3E(4)
      DOUBLE PRECISION PK1,PK2,PK3,PL
      COMMON /QQG/ PK1(4),PK2(4),PK3(4),PL(4)
      DOUBLE PRECISION QPR(4),PP(4)
      DOUBLE PRECISION a1ma1,Q2_c
      COMMON /QQG_CON/a1ma1,Q2_c
      DOUBLE PRECISION DBW2(4),DBMX(4)
      DOUBLE PRECISION STHETA,SPHI
      DOUBLE PRECISION MX2,M2,K12,K1K2,K22
      DOUBLE PRECISION WMAX
      INTEGER IMIX
      COMMON /OALPHAS/ WMAX,IMIX
      INTEGER IGENFL
      COMMON/GENWEI/IGENFL
      DOUBLE PRECISION XPOM
      COMMON/BARTELS/XPOM
      DOUBLE PRECISION WT1,sigdt_p,sigdt_m,sigdl,sigdi,WT2
      DOUBLE PRECISION alpha1,beta2,PREF,ALPH_EM,ALP,Y,W2,XB
      DOUBLE PRECISION m2bs,alpha2,beta1,a21ma2
      REAL SNGL
      REAL SLO,SHI,XEPS,EPS
      REAL SUMM11,XM11,SUMM21,XM21
      REAL SUMM12,XM12,SUMM22,XM22
      REAL SUMM1,XM1
      REAL SUMM2,XM2
      REAL FIL,FIU
      REAL L2MIN,L2MAX
      DOUBLE PRECISION DOT
      INTEGER I,iphase
      COMMON /QQG_C/L2MIN,L2MAX
      Double Precision   sst,betat,xppt,t,qqt
      common   /parameter/ sst,qqt,betat,xppt,t
      Double Precision DK,th,uh,qf2,k22fr
      EXTERNAL SUMM11,SUMM21,SUMM12,SUMM22,SUMM1,SUMM2
      EXTERNAL FIL,FIU,DOT,DK
      Double Precision D_XGX,glu,xglu,ksi,xcorr,alphas
      Double Precision FACQ2,beta,z,arg,Coeff,m2cc,mq2
	Double Precision stmax
	Common/uniglu/stmax

	Integer IGLU
      Common/gluon/Iglu
	
      EXTERNAL ALPHAS,D_XGX,XGLU
cc iphase = 1 only phase space
c      iphase = 1
c iphase = 0 x section
      iphase = 0
      t = dble(t2gki)
      WT1 = 0.D0
      WT2 = 0.d0
      XEPS = 0.05
      SLO = 0.
      SHI = 1.
c      L2MIN = 1.
      L2MIN = 0.3
	if(iglu.eq.5.or.iglu.eq.1) L2MIN = 0.00001
      L2MAX = 1000.
      Q2_c = Q2
      WT1 = 0.0D0
      y = DBLE(YY)
      Xb = Q2/Y/SSS
c         xb = Q2/(W2 + Q2)
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
      DBMX(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBMX(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBMX(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBMX(4)=  P(NIA1,4) + P(NIA1+1,4)
      MX2=DOT(DBMX,DBMX)
c  goto gamma p system
      DBW2(1)=  P(NIA1,1) + P(2,1)
      DBW2(2)=  P(NIA1,2) + P(2,2)
      DBW2(3)=  P(NIA1,3) + P(2,3)
      DBW2(4)=  P(NIA1,4) + P(2,4)
      W2=DOT(DBW2,DBW2)
c  goto gamma pom system
      DBW2(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBW2(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBW2(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBW2(4)=  P(NIA1,4) + P(NIA1+1,4)

      CALL PYROBO(0,0,0.D0,0.D0,
     +  -DBW2(1)/DBW2(4),-DBW2(2)/DBW2(4),
     +  -DBW2(3)/DBw2(4))
      SPHI = PYANGL(P(NIA1,1),P(NIA1,2))
      CALL PYROBO(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = PYANGL(P(NIA1,3),P(NIA1,1))
      CALL PYROBO(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
c      call pylist(1)

c new prefactor for the use with lorentz invariant phase space
c according to mark
      Pref = 9.d0/32.D0/PI
c sum over light charges = 2/3
      IF(IHFLA.LE.3) THEN
         QF2 = 2.D0/3.D0
      ELSEIF(IHFLA.EQ.4) THEN
         QF2 = 4.D0/9.D0
      ELSEIF(IHFLA.EQ.5) THEN
         QF2 = 1.D0/9.D0
      ENDIF
      Pref = Pref*QF2*ALPH_EM*ALP**3

      wt1 = alph_em/y/Q2/pi*( (1.d0 + (1.d0 - y)**2)/2.d0 *
     +sigdt_p - 2.d0*(1.d0 - y) * sigdt_m + (1.d0 - y) * sigdl + (2.d0
     +-y)*dsqrt(1.d0 - y)*sigdi)
      IF(WT1.LT.0) THEN
         wt1  = 0.
      ENDIF
      IF(WT1.ne.wt1) THEN
         write(6,*) ' error wt1  ',wt1
      ENDIF
      WT1 = 1d0
      CALL PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL PYROBO(0,0,0.D0,0.D0,
     +  DBW2(1)/DBW2(4),DBW2(2)/DBW2(4),
     +  DBW2(3)/DBW2(4))
      RETURN
      END
