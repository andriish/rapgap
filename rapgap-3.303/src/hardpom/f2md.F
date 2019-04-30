      SUBROUTINE F2MD(BETA,X_POM,Q2,T2,XDX,FLQ)
* calculate  F_2^D in the approach of M. Diehl
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
#include "rgparas.inc"
#include "rgparam.inc"
#include "rgluco.inc"
      DOUBLE PRECISION mq,mp,mx2,beta0,mu0
      DOUBLE PRECISION PT2GEN,PHIGEN
      DOUBLE PRECISION XDX(-4:4)
      COMMON/HARDPOM/PT2GEN,PHIGEN
      DOUBLE PRECISION PT2MD,PHIMD
      COMMON/PT/PT2MD,PHIMD
      COMMON /PARMD/ bet,xpom,q2t,t2t,mq,eq2
      COMMON /SEL/IPHI
      EXTERNAL PYCHGE
      REAL SLO,SHI,FMD,XEPS,RESULT
#include "rgludat2.inc"
      DATA SMALL/1.D-6/
      EXTERNAL FMD
      IPHI = 0
      pi=2.d0*dasin(1.d0)
      XEPS = 0.01
      SLO = SNGL(SMALL)
      SHI = 1. - SLO
      q2t = q2
      t2t = t2
      bet = beta
      xpom = x_pom
c      write(6,*) ' xpom',x_pom,beta,q2,t2
      mp = 0.938d0
      mx2 = Q2*(1.d0-beta)/beta + t2
c      write(6,*) ' f2md mx2,t2 ',mx2,t2
      eq2 = 2.d0/3.d0
      beta0 = 2.d0
      mu0 = 1.2d0
      mq = 0.3d0
      alphas0 = 1.d0
      NI = 25
      DO 10 I=-4,4
   10 XDX(I)=0.d0
      F2GF = 0.D0
      pt2min=pt2cut(12)
c      write(6,*) ' F2MD ptmin,pt2cut',pt2min,pt2cut(12)
c      write(6,*) 'F2MD: beta,xpom,q2,t2,f1t',beta,x_pom,q2,t2
      F2QMAX = 0.D0
      DO 20 IP=1,NFLAV
c current quark masses
c         mq=dble(pymass(ip))
c constituent quark masses
         mq=dble(PARF(100+IP))
         eq2 = dble(PYCHGE(IP))**2/9.0d0
c         write(6,*) ' mass ',mq,eq2,ip
         F2QTF =0.d0
         CALL INTGA(SLO,SHI,FMD,XEPS,RESULT)
c        write(6,*) ' after INTGA : ',RESULT
         F2QTF=DBLE(RESULT)
         XDX(IP)=F2QTF/eq2/2.D0
   20 CONTINUE
      XDX(-4)=XDX(4)
      XDX(-3)=XDX(3)
      XDX(-2)=XDX(2)
      XDX(-1)=XDX(1)
      XDX(0) = 0.d0
c      write(6,*) ' F2MD:F2QTF ',F2QTF,XDX
c      write(6,*) ' F2MD final : PT2GEN,v,beta,q2',pt2gen,v,beta,q2
      if(iwei.eq.1) then
         F2QMAX = 0.D0
         q2t = q2
         t2t = t2
         bet = beta
         xpom = x_pom
c      write(6,*) ' xpom',x_pom,beta,q2,t2
         mp = 0.938d0
         mx2 = Q2*(1.d0-beta)/beta + t2
         DO 30 I1 = 0,NI
            ip = KPA
c current quark masses
c           mq=dble(pymass(ip))
c constituent quark masses
            mq=dble(PARF(100+IP))
            eq2 = dble(PYCHGE(IP))**2/9.0d0
            rn2 = dfloat(I1)/dfloat(NI)
            if(rn2.eq.0.d0) rn2 = small
            if(rn2.eq.1.d0) rn2 = 1.d0 - small


            F2QT = dble(FMD(SNGL(rn2)))
c         FLQ =  sigl * wpt2
            F2TEST = F2QT
c         write(6,*) ' f2md 1st: F2QT FLQ ',F2QT,FLQ,rn2
c         write(6,*) ' pt2,wpt2,sigt',pt2,wpt2,sigt
            IF(F2TEST.GE.F2QMAX) F2QMAX = F2TEST
   30    CONTINUE
         F2QMAX = 2.0D0*F2QMAX

   40    CONTINUE
         IPHI = 1
         rn1 = draprn()
         rn2 = draprn()
c chi integration (0 - 2*pi)
         PHIMD = rn1*2.d0*pi
         if(rn2.gt.0.998d0) goto 40
         F2GEN = dble(FMD(SNGL(rn2)))
         PT2GEN = PT2MD
         PHIGEN = PHIMD
         F2RN = draprn()
c            write(6,*) ' F2MD: KPA = ',IP
c            write(6,*) ' F2MD : F2GEN F2MAX ',F2GEN,F2QMAX
c            write(6,*) ' f2MD: pt2,pl2,chi ',pt2,pl2,chi
         IF(F2GEN.GE.F2QMAX) THEN
            write(6,*) ' F2MD : F2GEN > F2MAX ',F2GEN,F2QMAX
         ELSEIF(F2GEN.LT.0.D0.OR.F2QMAX.LT.0.D0) THEN
            write(6,*) ' F2MD : F2GEN OR F2MAX < 0 ',F2GEN,F2QMAX
            write(6,*) ' F2MD : beta,x_pom,q2,t2',beta,x_pom,q2,t2

         ENDIF
c rotation with pt2gen is done in pyremn similar to primordial pt
         IF(F2QMAX*F2RN.GT.F2GEN) GOTO 40
c        write(6,*) ' F2MD final : PT2GEN,PHIGEN',pt2gen,phigen
      endif
      RETURN
      END

      FUNCTION FMD(X)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
#include "rgparas.inc"
      DOUBLE PRECISION mq,mx2
      DOUBLE PRECISION PT2MD,PHIMD
      COMMON/PT/PT2MD,PHIMD
      COMMON /PARMD/ beta,xpom,q2,t2,mq,eq2
      EXTERNAL PYCHGE
      REAL X,FMD
      pi=2.d0*dasin(1.d0)
      pt2min=pt2cut(12)
      mx2 = Q2*(1.d0-beta)/beta + t2
      FLQ = 0.D0
      F2Q = 0.D0
      FMD = 0.
c      write(6,*) 'F2MD: beta,xpom,q2,t2,f1t',beta,x_pom,q2,t2,f1t
      wchi = 2.d0*pi
c try 1/pt**4
      pt2max = mx2/4.d0 - 4.d0*mq**2
c         write(6,*) pt2min,pt2max
      if(pt2min.gt.pt2max) return
      rn2 = dble(x)
      pt2 = pt2max*pt2min/(pt2max + rn2*(pt2min-pt2max))
      pl2 = mx2/4.d0 - mq**2 - pt2
      wpt2 = pt2**2*(pt2max-pt2min)/pt2max/pt2min
c multply wpt2 with jacobian pt -> pl : pt/pl
c and account for divergence at high pt --> pl --> 0
      if(rn2.gt.0.998d0) return
      wpt2=wpt2*dsqrt(pt2)/dsqrt(pl2)
      call sigmd(pt2,sigt,sigl)
c now from sigma to F2, FL
      F2Q = sigt*Q2/4.d0/pi * wpt2
      FLQ = sigl*Q2/4.d0/pi * wpt2
      IF(F2Q.LE.0.0) THEN
         write(6,*) 'FMD<0 ',sigt,sigl,wpt2
      ENDIF
c      write(6,*) 'F2Q ',F2Q
      pt2md=pt2
      FMD = SNGL(F2Q)
      RETURN
      END
      SUBROUTINE SIGMD(pt2,sigt,sigl)
c note that sigt,sigl is here defined as q2/4/pi*sigt (sigl)
c is actually the structure funtion
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
#include "rgparas.inc"
#include "rgrapgki.inc"
      DOUBLE PRECISION lambda2,mq,mp,mx2,beta0,mu0
      DOUBLE PRECISION PT2MD,PHIMD
      COMMON/PT/PT2MD,PHIMD
      COMMON/HARDPOM/PT2GEN,PHIGEN
      COMMON /PARMD/ beta,xpom,q2,t2,mq,eq2
      COMMON /SEL/IPHI
      EXTERNAL PYCHGE
      sigt = 0.d0
      sigl = 0.d0
      mp = 0.938d0
      mp = 1.D0
      mx2 = Q2*(1.d0-beta)/beta + t2
      alpha_p=1.08d0 + t2/2.d0
      f1t = (4.d0*mp**2 - 2.8d0*t2)/(4.d0*mp**2 - t2)
      f1t = f1t/(1.d0 - t2/0.7d0)**2
c      write(6,*) ' F2MD t2: ',t2
      beta0 = 2.d0
      mu0 = 1.1d0
      betamu = beta0*mu0
      alphas0 = 1.d0
      pi=2.d0*dasin(1.d0)
c      write(6,*) 'sigMD: beta,xpom,q2,t2,f1t',beta,xpom,q2,t2,f1t
c      write(6,*) ' pt2 ',pt2
      pl2 = mx2/4.d0 - mq**2 - pt2
      if(pl2.le.0.0D0.or.pt2.le.0.0d0) return
      pl = dsqrt(pl2)
      pt = dsqrt(pt2)

c now F2MD
      lambda2 = (pt**2 + mq**2)/(1.d0 - beta)
      alpha_s=alphas(dsqrt(lambda2))
      if(alpha_s.lt.0.d0) alpha_s=0.d0
      if(alpha_s.ge.1.d0) alpha_s=1.d0
      xm = mx2/(pt2+mq**2)
      xp = pt2/(pt2+mq**2)
      aa = (1.d0-beta)*xp
      sigpp = xm**2*(1.d0 -2.d0/xm)*xp*(1.d0-aa)**2
      sigpp = sigpp+ 0.25d0*xm**2*mq**2/(pt2+mq**2)*(1.d0-2.d0* aa)**2
      sigl = Q2/mx2 *xm *(1.d0 -2.d0*aa)**2
      sigpm = 2.d0*xm*xp*(1.d0-aa)**2
      sigp0 = -2.d0*dsqrt(q2)/dsqrt(mx2)*xm*pl/pt*xp
      sigp0 = sigp0*(1.d0-aa)*(1.d0-2.d0*aa)
      sigt = sigpp
      comfac = 1.d0/dsqrt(mx2)*f1t**2 *xpom**(2.d0* (1.d0-
     +alpha_p))/xpom*27.d0/2.d0/pi**3 *eq2*alpha_s/alphas0*betamu**4/
     +(mx2+q2)**2 *(1.d0-beta)
c            write(6,*) 'SIGMD: pt2,comfac,sigt',pt2,comfac,sigt
      epsilon =(1.d0 - dble(yy))/(1.d0 - dble(yy) + dble(yy)**2/2.d0)
      if(iphi.eq.1) then
         sigt = comfac/2.D0/pi*(sigt + epsilon*sigl
     +    - epsilon*dcos(2.d0*phimd)*sigpm
     +    - dsqrt(2.d0*epsilon*(1.d0+epsilon))*dcos(phimd)*sigp0)
         sigl = comfac*sigl
      else
         sigt = comfac*(sigt+epsilon*sigl)
         sigl = comfac*sigl
      endif
      if(sigt.lt.0.0) then
         write(6,*) ' sigt<0.0 ',comfac,epsilon,phimd,sigpp
         write(6,*) ' sigt1 ',- epsilon*dcos(2.d0*phigen)*sigpm
         write(6,*) ' sigt2 ',- dsqrt(2.d0*epsilon*(1.d0+epsilon))
         write(6,*) ' sigt3 ',dcos(phimd)*sigp0
      endif
      RETURN
      END
