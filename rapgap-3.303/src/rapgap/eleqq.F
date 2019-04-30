      SUBROUTINE ELEQQ(WT1)
C
C    e p  ----> Q Q_BAR  p
C
C
      implicit none
#include "rgpart.inc"
#include "rglujets.inc"
#include "rgpara1.inc"
#include "rgpara.inc"
#include "rgluco.inc"
#include "rgrapgki.inc"
      DOUBLE PRECISION DBW2(4),DBMX(4)
      DOUBLE PRECISION STHETA,SPHI,sphi1
      DOUBLE PRECISION MX2,M2
      DOUBLE PRECISION WMAX
      DOUBLE PRECISION I_L,I_T
      INTEGER IMIX
      COMMON /OALPHAS/ WMAX,IMIX
      INTEGER IGENFL
      COMMON/GENWEI/IGENFL
      DOUBLE PRECISION XPOM
      COMMON/BARTELS/XPOM
      REAL SNGL
      DOUBLE PRECISION DOT
      INTEGER I
      DOUBLE PRECISION   sst,betat,xppt,t,qqt
      common   /parameter/ sst,qqt,betat,xppt,t
      DOUBLE PRECISION WT1,sigdt,sigdl,sigdi,sigda
      DOUBLE PRECISION ALPH_EM,Y,W2,XB,PT2
      DOUBLE PRECISION PTFACT,mf,ksi,ALPHA_S,ALPHAS
      DOUBLE PRECISION pt2max,charge,const
      DOUBLE PRECISION phi
      DOUBLE PRECISION cosp,cos2p,raphi

      common   /phi/       phi
      INTEGER IPHI
      COMMON /SEL/IPHI
      DOUBLE PRECISION D_XGX,glu,dglu,xglu

      DOUBLE PRECISION PK1,PK2,PL
      COMMON /QQ / PK1(4),PK2(4),PL(4)
      DOUBLE PRECISION xmx2,Q2_c,mf2
      COMMON /QQ_CON/xmx2,Q2_c,mf2
	Double Precision stmax
	Common/uniglu/stmax
      REAL L2MIN,L2MAX
      COMMON /QQG_C/L2MIN,L2MAX
	Integer IGLU
      Common/gluon/Iglu
      REAL SLO,SHI,XEPS,EPS
      REAL FIL,FIU
      REAL XIT,XIL,RES_T,RES_L

	Integer Iqq
	
      EXTERNAL ALPHAS,D_XGX,XGLU,raphi
      EXTERNAL FIL,FIU
      EXTERNAL XIT,XIL
	
	Iqq = 0
	
      t = dble(t2gki)
      sst = SSS
      qqt = q2
      xppt = xpom
      iphi=0
      XEPS = 0.05
      SLO = 0.
      SHI = 1.
c      L2MIN = 1.
      L2MIN = 0.3
	if(iglu.eq.5.or.iglu.eq.1) L2MIN = 0.00001
      L2MAX = 1000.
      Q2_c = Q2
      WT1 = 0.D0
      y = DBLE(YY)
      DBMX(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBMX(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBMX(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBMX(4)=  P(NIA1,4) + P(NIA1+1,4)
      MX2=DOT(DBMX,DBMX)
	xmx2=mx2
c  goto gamma p system
      DBW2(1)=  P(NIA1,1) + P(2,1)
      DBW2(2)=  P(NIA1,2) + P(2,2)
      DBW2(3)=  P(NIA1,3) + P(2,3)
      DBW2(4)=  P(NIA1,4) + P(2,4)
      W2=DOT(DBW2,DBW2)
      xb = Q2/(W2 + Q2)

      betat = xb/xpom

c  goto gamma pom system
      DBW2(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBW2(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBW2(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBW2(4)=  P(NIA1,4) + P(NIA1+1,4)
c     write(6,*) ' before boost '
c      call PYLIST(1)
      CALL PYROBO(0,0,0.D0,0.D0,
     +  -DBW2(1)/DBW2(4),-DBW2(2)/DBW2(4),
     +  -DBW2(3)/DBw2(4))
c     write(6,*) ' before phi '
c      call PYLIST(1)
      SPHI = PYANGL(P(NIA1,1),P(NIA1,2))
      CALL PYROBO(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
c     write(6,*) ' before theta '
c      call PYLIST(1)
      STHETA = PYANGL(P(NIA1,3),P(NIA1,1))
      CALL PYROBO(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
c     write(6,*) ' after boost '
c      call PYLIST(1)

      DO 10 I=1,4
         PK2(I) = P(NF1,I)
         PK1(I) = P(NF2,I)
   10 continue
      SPHI1 = PYANGL(P(1,1),P(1,2))
      call PYROBO(0,0,0.D0,-sphi1,0.d0,0.d0,0.d0)
      phi =  PYANGL(P(NF1,1),P(NF1,2))
c      phitt = raphi(1,nia2,nf1)
c     write(6,*) 'eleqqf  phitt,phit ',phittt,phitt,phit
      call PYROBO(0,0,0.D0,sphi1,0.d0,0.d0,0.d0)
      if(phi.ge.0.0) then
         PHITGKI = SNGL(PHI)
      else
         PHITGKI = SNGL(2.D0*PI+PHI)
      endif

      cosp = dcos(phi)
      cos2p = dcos(2.d0*phi)
      m2 = DOT(PK2,PK2) + 2.D0*DOT(PK2,PK1) + DOT(PK1,PK1)
      m2 = mx2
      PT2 = DBLE(PT2H)
      y = DBLE(YY)
      IF(p(nf1,5).lt.1.d0) then
         mf = 0.0d0
         charge = 2.d0/3.d0
      elseif(p(nf1,5).gt.1.d0.and.p(nf1,5).lt.4.d0) then
         mf = p(nf1,5)
         charge = 4.d0/9.d0
      endif
	mf2=mf**2
      ksi = (Q2 + mf**2/pt2*(Q2+m2))/m2
      if(iqq.eq.1) then
      glu = xglu(xpom,pt2*(1.d0+ksi))
      dglu = 0.d0
      dglu = d_xgx(xpom,pt2*(1.d0+ksi))
calculate I_T
      I_T = 2.d0/(1.d0 + ksi)**3 *(
     +      2.d0*ksi/pt2  * glu +
     +      (1.d0 -ksi -2.d0*ksi*dlog(ksi/(1.d0+ksi)))*
     +       dglu)
calculate I_L
      I_L = (pt2 + mf**2)*Q2/pt2/m2/(1.d0 + ksi)**3
     +      *((ksi - 1.d0)/pt2 * glu +
     +        (2.d0 +(1.d0-ksi)*dlog(ksi/(1.d0+ksi)))*
     +         dglu)
c      write(6,*) '   integerated ',I_T,I_L
      elseif(iqq.eq.0) then
c calculate integral over l
calculate I_T
         EPS = XEPS
         CALL GADAP2(SLO,SHI,FIL,FIU,XIT,EPS,RES_T)
	   I_T = -RES_T
calculate I_L
         EPS = XEPS
         CALL GADAP2(SLO,SHI,FIL,FIU,XIL,EPS,RES_L)
	   I_L = -RES_L
c      write(6,*) ' unintegerated ',I_T,I_L
	endif
      if((PT2+mf**2)/m2.GT.0.25D0) GOTO 20
c jacobian dM**2 --> d x_pom
c now in partdh included
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)
      Q2Q = pt2*(1.d0+ksi)
      ALPHA_S=ALPHAS(DSQRT(Q2Q))
      pt2max=m2/4.d0 - 4.d0*PYMASS(3)**2
      if(pt2.ge.pt2max) goto 20
      PTFACT = (PT2 + mf**2)/DSQRT(1.D0 - 4.d0*(PT2+mf**2)/m2)
      const = alph_em*pi**2*alpha_s**2*charge

      sigdt = 1.d0/m2**2/pt2/12.d0*const*
     +  PTFACT*( (1.d0 - 2.d0*(pt2+mf**2)/m2)*I_T**2
     + + mf**2*pt2*m2**2/(pt2+mf**2)/Q2**2 * I_L**2
     + )
      sigdl = 1.d0/m2**2/q2*4.d0/3.d0*ptfact*const*I_L**2
      sigdi = 1.d0/m2**2/dsqrt(pt2*q2)/3.d0*const*(PT2 + mf**2)*I_T*I_L
      sigda = 1.d0/m2**3/pt2/12.d0*ptfact*const*(PT2 + mf**2)*I_T**2
c     write(6,*) ' sigs ',sigdt,sigdl,sigdi,sigda
c     write(6,*) ' I_T,I_L ',I_T,I_L,ksi
c      wt1 = alph_em/2.d0/y/Q2/pi**2*(
      wt1 = alph_em/y/Q2/pi*(
     +      (1.d0 + (1.d0 - y)**2)/2.d0 * sigdt
     +       + (1.d0 - y) * sigdl - 2.d0*(1.d0-y)*cos2p*sigda
     +       + (2.d0 -y)*dsqrt(1.d0 - y)*cosp*sigdi)
c include t dependence a la DL
c now in partdh included
      IF(WT1.LT.0) THEN
         write(6,*) ' wt1  ',wt1
      ENDIF
      IF(WT1.NE.WT1) THEN
         write(6,*) ' wt1  ',wt1
         write(6,*) ' sigblw:sigdt,sigdl',sigdt,sigdl
         write(6,*) sigdt,sigdl,sigdi,sigda,ptfact,const,m2,I_T,I_L,
     +   PT2
         call PYLIST(1)
      ENDIF
      CALL PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL PYROBO(0,0,0.D0,0.D0,
     +  DBW2(1)/DBW2(4),DBW2(2)/DBW2(4),
     +  DBW2(3)/DBW2(4))
      RETURN
   20 WT1 = 0.D0
c      write(6,*) ' eleqq: set wt1=0'
      CALL PYROBO(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL PYROBO(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL PYROBO(0,0,0.D0,0.D0,
     +  DBW2(1)/DBW2(4),DBW2(2)/DBW2(4),
     +  DBW2(3)/DBW2(4))
      RETURN
      END
      FUNCTION XIT(X,PHI)
      implicit none
*KEEP,RGPARAM.
      INTEGER IWEI
      DOUBLE PRECISION ALPHS,PI,ALPH
      COMMON /RGPARAM/ ALPHS,PI,ALPH,IWEI
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
*KEND.
      REAL XIT,X,PHI
      DOUBLE PRECISION TIT,F_G,Y,WT
      DOUBLE PRECISION PK1,PK2,PL
      COMMON /QQ / PK1(4),PK2(4),PL(4)
      REAL L2MIN,L2MAX
      COMMON /QQG_C/L2MIN,L2MAX
      EXTERNAL TIL,F_G
      Y = DBLE(L2MIN*(L2MAX/L2MIN)**X)
      WT = Y *DBLE(LOG(L2MAX/L2MIN))
      WT = WT/2.D0/PI
c  phi dep
      PL(1) = DSQRT(Y)*DBLE(COS(PHI))
      PL(2) = DSQRT(Y)*DBLE(SIN(PHI))
      XIT = sngl(WT*TIT(1)*F_G(Y)/Y)
c          write(6,*) ' xit x,phi ',xit,x,phi
c      IF(XIT.LT.0.) THEN
c         write(6,*) ' xit x,phi ',x,phi
c         write(6,*) 'TIT(1,1),F_G(Y) ',TIT(1),F_G(Y)
c      ENDIF
      RETURN
      END
      FUNCTION XIL(X,PHI)
      implicit none
#include "rgparam.inc"
      REAL XIL,X,PHI
      DOUBLE PRECISION TIT,F_G,Y,WT
      DOUBLE PRECISION PK1,PK2,PL
      COMMON /QQ / PK1(4),PK2(4),PL(4)
      REAL L2MIN,L2MAX
      COMMON /QQG_C/L2MIN,L2MAX
      EXTERNAL TIL,F_G
      Y = DBLE(L2MIN*(L2MAX/L2MIN)**X)
      WT = Y *DBLE(LOG(L2MAX/L2MIN))
      WT = WT/2.D0/PI
c  phi dep
      PL(1) = DSQRT(Y)*DBLE(COS(PHI))
      PL(2) = DSQRT(Y)*DBLE(SIN(PHI))
      XIL = sngl(WT*TIT(2)*F_G(Y)/Y)
c          write(6,*) ' summ11 x,phi ',summ11,x,phi
      IF(XIL.LT.0.) THEN
c         write(6,*) ' summ11 x,phi ',x,phi
c         write(6,*) 'TIL(1,1),F_G(Y) ',TIL(1,1),F_G(Y)
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION TIT(I)
      implicit none
      DOUBLE PRECISION PK1,PK2,PL
      COMMON /QQ / PK1(4),PK2(4),PL(4)
      INTEGER I
      DOUBLE PRECISION xmx2,Q2_c,mf2
      COMMON /QQ_CON/xmx2,Q2_c,mf2
	
	Double Precision TILT,XK2,XL2,PKL
	
      TILT = 0.D0
	XK2 = PK1(1)**2 + PK1(2)**2
	XL2 = PL(1)**2 + PL(2)**2
	PKL = (xl2+xk2/xmx2*(Q2_c-xmx2) + mf2/xmx2*(Q2_c+xmx2))**2
     +      + 4d0*xk2*(xk2*Q2_c/xmx2 + mf2/xmx2*(Q2_c+xmx2))
	IF(I.EQ.1) THEN
c transverse part (I_T)	
c	TILT = (xmx2 - Q2_c)/(xmx2+Q2_c) +
c     +(xl2 + xk2/xmx2*(Q2_c-xmx2))/dsqrt(PKL)
	TILT = 2.*xmx2*xk2 /(xmx2+Q2_c)/(xk2+mf2)-1. +
     +(xl2 + xk2/xmx2*(Q2_c-xmx2)+mf2*(1.+Q2_c/xmx2))/dsqrt(PKL)
	ELSEIF(I.EQ.2) THEN
c longitudinal part (I_L)	
	TILT = Q2_c/(xmx2+Q2_c) -
     +	 (xk2+mf2)*Q2_c/xmx2/dsqrt(PKL)
	
	ELSE
	write(6,*) ' TIT: I not possible, I = ',I
	STOP
	Endif
c	write(6,*) ' TIT: Q2_c,xmx2,xl2,xk2 ',Q2_c,xmx2,xl2,xk2

c      IF(TILT.LT.0) THEN
c      write(6,*) ' TIT I Q2_c,xmx2,xl2,xk2 ',TILT,I,Q2_c,xmx2,xl2,xk2
c      ENDIF
      TIT = TILT
      RETURN
      END
