      SUBROUTINE F2BLW(BET,X_POM,Q2,T2,F2QT,FLQ)
* calculate  F_2^D in the hard approach of J. Bartels, H. Lotter, M. Wuesthoff
*                                          (DESY 96-026)
* the subroutines come from H. Lotter
* THANKS A LOT
*

      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
      common   /phi/       phi
      common   /parameter/ s,q,beta,xp,t
      INTEGER IPHI
      COMMON /SEL/IPHI
#include "rgpart.inc"
#include "rgpara.inc"

      REAL * 4 SLO,SHI,FBLW,XEPS,RESULT,PYMASS
      EXTERNAL FBLW,PYMASS
      DATA SMALL/1.D-6/
      IPHI = 0
      XEPS = 0.01
      SLO = SNGL(SMALL)
      SHI = 1. - SLO
      q = Q2
      beta = bet
      xp = x_pom
      t = t2
      phi = 0.d0
      aem=1.d0/137.d0
      pi=2.d0*dasin(1.d0)
      s = SSS
      pt2min=pt2cut(12)
      xm=q2/beta-q2 + t
      pt2max=xm/4.d0 - 4.d0*DBLE(PYMASS(3))**2
      F2QT = 0.D0
      FLQ = 0.D0
      if(pt2min.gt.pt2max) return
      NI = 25

c      write(6,*) 'F2BLW: beta,xpom,q2,t2',beta,x_pom,q2,t2
      CALL INTGA(SLO,SHI,FBLW,XEPS,RESULT)
c      write(6,*) ' after INTGA : ',RESULT
c          write(6,*) ' f2blw  xng,xpom',,xng,xpom
c      write(6,*) 'f2blw f2qt ',f2qt
      F2QT=DBLE(RESULT)
      FLQ=0.0D0

   10 CONTINUE
c         write(6,*) ' f2blw iwei =',iwei
      if(iwei.eq.1) then
c         write(6,*) ' f2blw iwei =1'
         F2QMAX = 0.D0
         DO 20 I1 = 0,NI
            rn1 = dfloat(I1)/dfloat(NI)
            if(rn1.eq.0.d0) rn1 = small
            if(rn1.eq.1.d0) rn1 = 1.d0 - small
c weighting with 1/pt**4
            pt2 = pt2min*pt2max/(rn1*pt2min + pt2max*(1.D0 - rn1))
            wpt2 = pt2**2*(pt2max-pt2min)/pt2min/pt2max
            if(rn1.gt.0.998d0) goto 20
            CALL SIGBLW(pt2,sigt,sigl)
            F2TEST = sigt*Q2/(4.D0*PI**2*aem)*wpt2
c         write(6,*) ' f2blw 1st: F2TEST,wpt2 ',F2TEST,wpt2,sigt,pt2
            IF(F2TEST.GE.F2QMAX) F2QMAX = F2TEST

   20    CONTINUE
         F2QMAX = 2.D0*F2QMAX
   30    CONTINUE
         iphi = 1
         rn1 = draprn()
c weighting with 1/pt**4
         pt2 = pt2min*pt2max/(rn1*pt2min + pt2max*(1.D0 - rn1))
         wpt2 =  pt2**2*(pt2max-pt2min)/pt2min/pt2max
         if(rn1.gt.0.998d0) goto 30
         rn2 = draprn()
c chi integration (0 - 2*pi)
         PHI = rn2*2.d0*pi

         CALL SIGBLW(pt2,sigt,sigl)
         F2GEN =  sigt*Q2/(4.D0*PI**2*aem)*wpt2
         F2RN = draprn()
         IF(F2GEN.GE.F2QMAX) THEN
            write(6,*) ' F2BLW : F2GEN > F2MAX ',F2GEN,F2QMAX
         ELSEIF(F2GEN.LT.0.D0.OR.F2QMAX.LT.0.D0) THEN
            write(6,*) ' F2BLW: F2GEN OR F2MAX < 0 ',F2GEN,F2QMAX
            write(6,*) ' F2BLW: beta,x_pom,q2,t2',beta,x_pom,q2,t2
            write(6,*) ' F2BLW: sigt,wpt2 ',sigt,wpt2
         ENDIF
c rotation with pt2gen is done in pyremn similar to primordial pt
         IF(F2QMAX*F2RN.GT.F2GEN) THEN
c            write(6,*) 'F2QMAX*F2RN.GT.F2GEN',F2QMAX,F2RN,F2GEN
            GOTO 30
         ENDIF
         PT2GEN = pt2
         PHIGEN = PHI
c         write(6,*) ' F2BLW final : PT2GEN',pt2gen,' phigen ',phigen
      endif
      RETURN
      END

      FUNCTION FBLW(X)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
      common   /phi/       phi
      common   /parameter/ s,q,beta,xp,t
#include "rgpart.inc"
#include "rgparas.inc"
      REAL * 4 X,FBLW,PYMASS
      EXTERNAL PYMASS
      aem=1.d0/137.d0
      pi=2.d0*dasin(1.d0)
      pt2min=pt2cut(12)
      xm=q/beta-q + t
      pt2max=xm/4.d0 - 4.d0*DBLE(PYMASS(3))**2
      F2QT = 0.D0

      FBLW = 0.0
c      write(6,*) ' FBLW ',pt2min,pt2max,x
      if(pt2min.gt.pt2max) return
      rn1 = dble(x)
c weighting with 1/pt**2
c      pt2 = pt2min*(pt2max/pt2min)**rn1
c      wpt2 = log(pt2max/pt2min)*pt2
c weighting with 1/pt**4
      pt2 = pt2max*pt2min/(pt2max + rn1*(pt2min-pt2max))
      wpt2 = pt2**2*(pt2max-pt2min)/pt2max/pt2min
      CALL SIGBLW(pt2,sigt,sigl)
c         write(6,*) 'f2blw sig,pt2 ',f2qt,pt2
      F2QT=sigt*q/(4.D0*PI**2*aem)*wpt2
      FLQ=sigL*q/(4.D0*PI**2*aem)*wpt2
c          write(6,*) ' f2blw xng,xpom',,xng,xpom
c       write(6,*) ' s,q,beta,xp,t ',s,q,beta,xp,t
c       write(6,*) 'f2blw f2qt ',f2qt,pt2,x
      FBLW = SNGL(F2QT)
      RETURN
      END

      subroutine sigblw(pt,sigt,sigl)
c
c     <sig> berechnet den wq als funktion von pt^2, die anderen
c     kinematischen variablen werden im common uebergeben.
c     pt,q und m stehen immer fuer physikalisch pt^2,Q^2,M^2 .
c
c     die gluonstrukturfunktion wird im up grv_nl berechnet,
c     die ableitung der strfkt. in dgrv_nl.
c
c     das up gamma berechnet einen vorfaktor. es ruft ein
c     programm zur berechnung von running alpha_s auf.
c     durch das up dlformf wird implizit die t-integration
c     ausgefuehrt.
c
c     sigt, sigl, und sigi repraesentieren den gamma^* - p
c     wq fuer transversale und longitudinale photonen sowie
c     den interferenzterm.
c
      implicit none
      DOUBLE PRECISION   sigt,sigl,sigi,sigtt,epsilon
      DOUBLE PRECISION   dlformf
      DOUBLE PRECISION   s,q,xb,pt,m,beta,xp,t
      DOUBLE PRECISION   fact,facl,pi
      DOUBLE PRECISION   grv_nl,dgrv_nl
      DOUBLE PRECISION   rggamma
      DOUBLE PRECISION   phi
c
      common   /phi/       phi
      common   /parameter/ s,q,beta,xp,t
      INTEGER IPHI
      COMMON /SEL/IPHI
#include "rgrapgki.inc"
c
      pi=2.d0*dasin(1.d0)
      fact=4.d0/3.d0
      facl=4.d0/3.d0
c
      m=q/beta-q
      xb=beta*xp
c
c
      sigtt=rggamma(pt*(q+m)/m)*
     +      fact*
     +     (1.d0-2.d0/m*pt)/dsqrt(1.d0-4.d0*pt/m)/pt**2*
     +     m**2/(q+m)**4*
     +     (
     +      grv_nl( xb*(q+m)/q,pt*(q+m)/m)*beta    +
     +      1.d0*
     +      pt*dgrv_nl( xb*(q+m)/q,pt*(q+m)/m)*(
     +       1.5d0*(1.d0-4.d0/3.d0*beta)  -
     +       1.d0*(1.d0-beta-(1.d0-beta)*dlog(beta)+dlog(beta))  )
     +     )**2
     +     *dlformf(m,q,xb,0)
c include jacobian dM**2 --> dx_pom
      sigtt = sigtt*q/xb
c
c
      sigl=rggamma(pt*(q+m)/m)*
     +     facl*1.d0/dsqrt(1.d0-4.d0*pt/m)*
     +     q/(q+m)**4/pt*
     +     (
     +      grv_nl( xb*(q+m)/q,pt*(q+m)/m)*(2.d0*beta-1.d0) +
     +      1.d0*
     +      pt*dgrv_nl( xb*(q+m)/q,pt*(q+m)/m)*(
     +       3.d0*(1.d0-beta)  -
     +       (1.d0-beta-2.d0*(1.d0-beta)*dlog(beta)+dlog(beta))  )
     +     )**2
     +     *dlformf(m,q,xb,0)
c include jacobian dM**2 --> dbeta
      sigl = sigl*q/xb
c
      sigi=rggamma(pt*(q+m)/m)
     +      /m**2*dsqrt(pt/q)*4.d0/3.d0*
     +      m*q*m**2/pt**2/(q+m)**4*
     +     (
     +      grv_nl( xb*(q+m)/q,pt*(q+m)/m)*beta    +
     +      1.d0*
     +      pt*dgrv_nl( xb*(q+m)/q,pt*(q+m)/m)*(
     +       1.5d0*(1.d0-4.d0/3.d0*beta)  -
     +       1.d0*(1.d0-beta-(1.d0-beta)*dlog(beta)+dlog(beta))  )
     +     )*
     +     (
     +      grv_nl( xb*(q+m)/q,pt*(q+m)/m)*(2.d0*beta-1.d0) +
     +      1.d0*
     +      pt*dgrv_nl( xb*(q+m)/q,pt*(q+m)/m)*(
     +       3.d0*(1.d0-beta)  -
     +       (1.d0-beta-2.d0*(1.d0-beta)*dlog(beta)+dlog(beta))  )
     +     )
     +     *dlformf(m,q,xb,0)
c include jacobian dM**2 --> dbeta
      sigi = sigi*q/xb
      epsilon =(1.d0 - dble(yy))/(1.d0 - dble(yy) + dble(yy)**2/2.d0)
      if(iphi.eq.1) then
         sigt =    ((1.d0-2.d0*epsilon*dcos(2.d0*phi)*1.d0*
     +             pt/m/(1.d0-2.d0*pt/m)
     +              )*sigtt
     +             + (2.d0-dble(yy))*dsqrt(1.d0-dble(yy))*sigi*dcos(phi)
     +               /(1.d0 - dble(yy) + dble(yy)**2/2.d0)
     +             + epsilon*sigl
     +              )
     +              /2.d0/pi
         sigl = sigl/2.d0/pi
      else
         sigt = sigtt + epsilon*sigl
      endif
      if(sigt.lt.0.0) then
         write(6,*) ' sigt<0.0',sigtt,epsilon
         write(6,*) ' sigt1',-2.d0*epsilon*dcos(2.d0*phi)*1.d0*
     +    pt/m/(1.d0-2.d0*pt/m)
         write(6,*) ' sigt2',(2.d0-dble(yy))*dsqrt(1.d0-dble(yy))
     +    *sigi*dcos(phi)
         write(6,*) ' sigt3 yy',yy,q/xb/s,xb,phi
         write(6,*) ' sigl sigt',sigl,sigt
      endif
c
      return
      end


      function rggamma(q)
      DOUBLE PRECISION   rggamma,q
      DOUBLE PRECISION   aem,alphasl,pi,charge
      common   /cons/   aem,pi
comment on aem: since we use f2 aem is cancelled, so it doesnt matter
      aem=1.d0/137.d0
      pi=2.d0*dasin(1.d0)
      charge=1.D0/3.D0*1.D0/3.D0+2.D0/3.D0*2.D0/3.D0+1.D0/3.D0*1.D0/3.D0
      rggamma=aem*alphasl(q)**2*pi**2*charge
      return
      end
c
      FUNCTION ALPHAsl(Q)
      IMPLICIT REAL*8 (A-H,P-Z)
      pi=2.d0*dasin(1.d0)
      ZETAM = 91.161D0
      ALPHAZ = 0.12D0
      FLAVOR = 4.0D0
      BZERO = 11.0D0 - 2.0D0*FLAVOR/ 3.0D0
      ALPHAsl = ALPHAZ/ (1 + ALPHAZ*BZERO*DLOG(dsqrt(Q)/ZETAM)/2.0D0/PI)
c here use for consistency alphas(SQRT(Q))
c      write(6,*) ' blw alphas ',alphasl
      alphasl = alphas(dsqrt(Q))
c      write(6,*) ' my alphas ',alphasl
      RETURN
      END


      function dlformf(m,q,xb,i)
      DOUBLE PRECISION   dlformf,m,xb,q,xp
      DOUBLE PRECISION   a,b1,b2
      DOUBLE PRECISION   dexpin
      DOUBLE PRECISION   s,wmax,wmin
      DOUBLE PRECISION   ss,beta,xpp,t,qq,alphap,dl
      common   /parameter/ ss,qq,beta,xpp,t
      integer  i
c
      common   /hera/   s,wmax,wmin
c
      xp=xb*(q+m)/q
      a=dlog(1.d0/xp)
c
      if(i.eq.1) then
         dlformf=0.30119d0
      else if(i.eq.2) then
         b1=dexpin(0.35d0*a)
         b2=dexp(0.35d0*a)
         dlformf=-0.00191*b1*b2*a**3-0.085699*a+0.04559*dexp(2.d0*a)*
     +   dexpin(2.d0*a)+0.00546*a**2-0.052477*a*dexp(2*a)* dexpin(2.d0*
     +   a)-0.02274*a*b1*b2+0.0245331*a**2*b1*b2+ 0.380649-0.04558633*
     +   b1*b2
      elseif(i.eq.0) then
c use donnachie landshoff formfactor
c        f(t)=(1/x)^(alpha'*t)*(4-2.8*t)/(4-t)*1/(1-t/0.7)^2
c        alpha' = 0.25
         alphap = 0.25d0
         dl =(1.d0/xp)**(alphap*t)*(4.d0-2.8d0*t)/(4.d0-t)*
     +   1.d0/(1.d0-t/0.7d0)**2
c         write(6,*) ' dlformf : ',dl
         dlformf = dl**2
      end if
      return
      end


      function grv_nl(x,q)
      implicit none
      DOUBLE PRECISION grv_nl,x,q,s1,s2
      s1 = (1-x)**(0.3404D1+0.83D0*dlog(0.5848488801D0*dlog(0.162591051D
     +2*q)))
      s2 = x**(0.1724D1+0.157D0*dlog(0.5848488801D0*dlog(0.162591051D2*q
     +)))*(0.7517D1-0.2547D1*dlog(0.5848488801D0*dlog(0.162591051D2*q))+
     +(0.3409D2-0.5221D2*dsqrt(dlog(0.5848488801D0*dlog(0.162591051D2*q)
     +))+0.1747D2*dlog(0.5848488801D0*dlog(0.162591051D2*q)))*x+(0.4039D
     +1+0.1491D1*dlog(0.5848488801D0*dlog(0.162591051D2*q)))*x**2)*dlog(
     +1/x)**(0.8D0+0.1016D1*dlog(0.5848488801D0*dlog(0.162591051D2*q)))+
     +dlog(0.5848488801D0*dlog(0.162591051D2*q))**0.1014D1*dexp(0.1112D1
     +-0.3438D1*dlog(0.5848488801D0*dlog(0.162591051D2*q))+0.302D0*dlog(
     +0.5848488801D0*dlog(0.162591051D2*q))**2+dsqrt((0.3256D1-0.436D0*d
     +log(0.5848488801D0*dlog(0.162591051D2*q)))*dlog(0.5848488801D0*dlo
     +g(0.162591051D2*q))**0.1738D1*dlog(1/x)))
      grv_nl = s1*s2
      return
      end


      function dgrv_nl(x,q)
      implicit none
      DOUBLE PRECISION   dgrv_nl,x,q
      DOUBLE PRECISION s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      DOUBLE PRECISION s11,s12,s13,s14
      integer itest
      data itest/1/
      if(itest.eq.0) then
         dgrv_nl = 0.d0
      else
         s2 = 0.83D0*(1-x)**(0.3404D1+0.83D0*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q)))/q
         s3 = 1/dlog(0.162591051D2*q)*dlog(1-x)*(x**(0.1724D1+0.157D0*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q)))*(0.7517D1-
     +   0.2547D1*dlog(0.5848488801D0*dlog(0.162591051D2*q))+(0.3409D2-
     +   0.5221D2*dsqrt(dlog(0.5848488801D0*dlog(0.162591051D2*q)))+
     +   0.1747D2*dlog(0.5848488801D0*dlog(0.162591051D2*q)))*x+
     +   (0.4039D1+0.1491D1*dlog(0.5848488801D0*dlog(0.162591051D2*q)))
     +   *x**2)*dlog(1/x)**(0.8D0+0.1016D1*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q)))+dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q))**0.1014D1*dexp(0.1112D1-0.3438D1*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q))+0.302D0*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q))**2+dsqrt((0.3256D1-
     +   0.436D0*dlog(0.5848488801D0*dlog(0.162591051D2*q)))*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q))**0.1738D1*dlog(1/x)
     +   )))
         s1 = s2*s3
         s3 = (1-x)**(0.3404D1+0.83D0*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q)))
         s6 = 0.157D0*x**(0.1724D1+0.157D0*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q)))/q/dlog(0.162591051D2*q)*dlog(x)*
     +   (0.7517D1-0.2547D1*dlog(0.5848488801D0*dlog(0.162591051D2*q))+
     +   (0.3409D2-0.5221D2*dsqrt(dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q)))+0.1747D2*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q)))*x+(0.4039D1+0.1491D1*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q)))*x**2)*dlog(1/x)**
     +   (0.8D0+0.1016D1*dlog(0.5848488801D0*dlog(0.162591051D2*q)))
         s7 = x**(0.1724D1+0.157D0*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q)))*(-0.2547D1/q/dlog(0.162591051D2*q)+(-
     +   0.26105D2/dsqrt(dlog(0.5848488801D0*dlog(0.162591051D2*q)))/q/
     +   dlog(0.162591051D2*q)+0.1747D2/q/dlog(0.162591051D2*q))*x+
     +   0.1491D1/q/dlog(0.162591051D2*q)*x**2)*dlog(1/x)**(0.8D0+
     +   0.1016D1*dlog(0.5848488801D0*dlog(0.162591051D2*q)))
         s5 = s6+s7
         s6 = s5+0.1016D1*x**(0.1724D1+0.157D0*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q)))*(0.7517D1-0.2547D1*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q))+(0.3409D2-0.5221D2*
     +   dsqrt(dlog(0.5848488801D0*dlog(0.162591051D2*q)))+0.1747D2*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q)))*x+(0.4039D1+
     +   0.1491D1*dlog(0.5848488801D0*dlog(0.162591051D2*q)))*x**2)*
     +   dlog(1/x)**(0.8D0+0.1016D1*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q)))/q/dlog(0.162591051D2*q)*dlog(dlog(1/x)
     +   )
         s7 = s6
         s9 = 0.1014D1*dlog(0.5848488801D0*dlog(0.162591051D2*q))**
     +   0.14D-1*dexp(0.1112D1-0.3438D1*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q))+0.302D0*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q))**2+dsqrt((0.3256D1-0.436D0*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q)))*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q))**0.1738D1*dlog(1/x)
     +   ))/q/dlog(0.162591051D2*q)
         s11 = dlog(0.5848488801D0*dlog(0.162591051D2*q))**0.1014D1
         s13 = -0.3438D1/q/dlog(0.162591051D2*q)+0.604D0*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q))/q/
     +   dlog(0.162591051D2*q)+1/dsqrt((0.3256D1-0.436D0*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q)))*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q))**0.1738D1*dlog(1/x)
     +   )*(-0.436D0/q/dlog(0.162591051D2*q)*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q))**0.1738D1*dlog(1/x)+0.1738D1*(0.3256D1-
     +   0.436D0*dlog(0.5848488801D0*dlog(0.162591051D2*q)))*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q))**0.738D0*dlog(1/x)/
     +   q/dlog(0.162591051D2*q))/2
         s14 = dexp(0.1112D1-0.3438D1*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q))+0.302D0*dlog(0.5848488801D0*
     +   dlog(0.162591051D2*q))**2+dsqrt((0.3256D1-0.436D0*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q)))*
     +   dlog(0.5848488801D0*dlog(0.162591051D2*q))**0.1738D1*dlog(1/x)
     +   ))
         s12 = s13*s14
         s10 = s11*s12
         s8 = s9+s10
         s4 = s7+s8
         s2 = s3*s4
c multiplying with q give p_t**2 *dgrv	
c         dgrv_nl = (s1+s2)*q
c now put this factor in x section calc...
         dgrv_nl = (s1+s2)
      endif
      return
      end
