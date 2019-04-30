      SUBROUTINE F2DHW(BETA,X_POM,Q2,T2,F2QTF,FLQF,xgl)
* calculate  F_2^D in the hard approach of M. Wuesthoff
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
#include "rgpara.inc"
      EXTERNAL POM,FG
      DATA SMALL/1.D-6/
      q02=1.0d0
      NTRY = 100
      NI = 5
      F2QTF = 0.D0
      FLQF = 0.D0
      xgl = 0.D0
c      write(6,*) 'F2DHW: beta,xpom,q2,t2',beta,x_pom,q2,t2
c      write(6,*) ' f2dhw iwei',iwei
      xi0=log(4.0d0*q02/q2*beta)
      xi1=12.d0
      F2QMAX = 0.D0
      DO 10 IRN = 1,NTRY
         rn1 = draprn()
         v=exp(xi0+rn1*(xi1-xi0))
         wv=(xi1-xi0)*v
c
         rn2 = draprn()
         z=beta/(beta+rn2*(1.d0-beta))
         wz = z**2 *(1.d0 -beta)/beta
         CALL XNTOTAL(F2QT,F2G,F2QL,Q2,x_pom,BETA,V,Z,T2)
         F2QTF = F2QTF + F2QT * wv/NTRY + F2G * wv * wz/NTRY
         FLQF =FLQF + F2QL * wv/NTRY
c get the gluon
         xng = 0.d0
c
c
         xpom=pom(x_pom,-t2,v)
         xpom = xpom/x_pom
         if (v.lt.q02/q2) then
            xng=0.0d0
         else
            xng=fg(beta,v)/v
            xng=xng**2/q02
            xng=xng*9.0d0/4.0d0
c
            xng = xng/16.d0 * xpom
         endif
         xgl = xgl + xng*wv/ntry

   10 CONTINUE
c          write(6,*) ' f2dhw f2qtf xng,xpom',f2qtf,xng,xpom

c      write(6,*) 'f2dhw : ',f2gf
      DO 20 I1 = 0,NI
         DO 20 I2 = 0,NI
            rn1 = dfloat(I1)/dfloat(NI)
            if(rn1.eq.0.d0) rn1 = small
            if(rn1.eq.1.d0) rn1 = 1.d0 - small
            v=exp(xi0+rn1*(xi1-xi0))
            wv=(xi1-xi0)*v
c
            rn2 = dfloat(I2)/dfloat(NI)
            if(rn2.eq.0.d0) rn2 = small
            if(rn2.eq.1.d0) rn2 = 1.d0 - small
            z=beta/(beta+rn2*(1.d0-beta))
            wz = z**2 *(1.d0 -beta)/beta
            CALL XNTOTAL(F2QT,F2G,F2QL,Q2,x_pom,BETA,V,Z,T2)
            F2TEST = (F2QT+F2QL)*wv + F2G*wv*wz
c         write(6,*) ' f2dhw 1st: (F2g) ',F2g,z,v,wv,wz
c         write(6,*) ' f2dhw 1st: (F2QT+F2QL),wv ',F2TEST,(F2QT+F2QL)*wv
            IF(F2TEST.GE.F2QMAX) F2QMAX = F2TEST

   20 CONTINUE
      F2QMAX = 5.D0*F2QMAX
   30 CONTINUE
      if(iwei.eq.1) then
         rn1 = draprn()
         v=exp(xi0+rn1*(xi1-xi0))
         wv=(xi1-xi0)*v
         rn2 = draprn()
         z=beta/(beta+rn2*(1.d0-beta))
         wz = z**2 *(1.d0 -beta)/beta
         CALL XNTOTAL(F2QT,F2G,F2QL,Q2,x_pom,BETA,V,Z,T2)
         F2GEN = (F2QT+F2QL)*wv + F2G*wv*wz
         F2RN = draprn()
         IF(F2GEN.GE.F2QMAX) THEN
            write(6,*) ' F2DHW : F2GEN > F2MAX ',F2GEN,F2QMAX
         ELSEIF(F2GEN.LE.0.D0.OR.F2QMAX.LE.0.D0) THEN
            write(6,*) ' F2DHW: F2GEN OR F2MAX < 0 ',F2GEN,F2QMAX
            write(6,*) ' F2DHW: beta,x_pom,q2,t2',beta,x_pom,q2,t2
         ENDIF
c rotation with pt2gen is done in pyremn similar to primordial pt
         IF(F2QMAX*F2RN.GT.F2GEN) GOTO 30
         FT = F2G*wv*wz/F2GEN
         PT2TT = (1.d0 - beta)*q02/v
         IF(FT.GE.draprn()) THEN
            PT2GEN = (1.d0 - z)*q02/v
         ELSE
            PT2GEN = (1.d0 - beta)*q02/v
         ENDIF
c      write(6,*) ' F2DHW final : PT2GEN,v,beta,q2,z',pt2gen,v,beta,q2,z
c      write(6,*) ' F2DHW final : PT2GEN(beta),PT2GEN',pt2tt,pt2gen
c      write(6,*) ' F2DHW final :F2GEN,F2G*wv*wz',F2GEN,F2G*wv*wz
      endif
      RETURN
      END
c-----------------------------------------------------
      double precision function pom(xp,t,v)
      implicit DOUBLE PRECISION (a-g,o-z)

c
c
      q02=1.0d0
      pm=1.09d0
      eps=0.08d0
      if (log(q02/v).gt.-2.0d0) then
         pm=pm/(log(q02/v)+3.0d0)**0.558d0
         eps=eps+0.0997d0*log(log(q02/v)+3.0d0)
      endif
      pm=pm**2
      eps=2.0d0*eps
      pm=pm*(0.05d0/xp)**eps
c
      fm=(4.0d0+2.8d0*t)**2/(4.0d0+t)**2
      form=fm/(1.0d0+t/0.7d0)**4*xp**(0.5d0*t)
c      form=exp(-2.0d0*b*t)*xp**(0.5d0*t)
c
      pom=pm*form
c      write(6,*) 'pom: xp,t,v ',xp,t,v
c      write(6,*) 'pom: pm,form',pm,form
c      write(6,*) ' old pom = ',pm*form
C new parameters from ANL-HEP-PR 97-03
c new pomeron flux (gluon density)
      A = 0.877d0
      B = 0.133d0
      C = 0.596d0
      x0 = 0.05d0
      alpom = 0.085
      scale = 1.d0/v
      lnscale=log(scale/q02)
      if(lnscale.lt.0d0) lnscale=0.d0
      if(scale.gt.q02) then
         alpom = alpom + B*log(lnscale+1.d0)
      endif
c here we use (x0/xp)**alpom instead of (x0/xp)**(1-alpom) as in the paper
c which is presumably a typing error
      pm = A * (x0/xp)**(alpom) *(lnscale + 1.d0)**(-C)
c new form of t dependence
      B_D = 5.9d0
      form = exp(-B_D*t)
c      form = form/B_D
c      write(6,*) ' wuesthoff v,t,xp:',v,t,xp,pm,form
      pom = pm**2 * form
c      write(6,*) ' new pom = ',pm*form

c
      return
      end
c------------------------------------------------------------
      double precision function fql(beta,v)
      implicit DOUBLE PRECISION (a-g,o-z)
      DOUBLE PRECISION nql,lnln

c
c      double precision beta,v,root,lnln,nql
c
      if (v.gt.1.0d3) then
         nql=log(beta*v)
      else
         root=sqrt(v**2+2.0d0*(1.0d0-2.0d0*beta)*v+1.0d0)
c
         lnln=root-(1.0d0-2.0d0*beta+v)
         lnln=lnln/(root+(1.0d0-2.0d0*beta)*v+1.0d0)
         lnln=log(lnln)
c
         nql=log(1.0d0/beta)
         nql=nql+lnln/root
         nql=nql+(1.0d0/root-1.0d0)*log(v)
      endif
c
      fql=nql
c       write(6,*)'nql=',nql,'v = ',v,' beta = ',beta
      return
      end
c---------------------------------------------------------------
      double precision function fqt(beta,v)
      implicit DOUBLE PRECISION (a-g,o-z)
      DOUBLE PRECISION nqt,lnln
c
c      double precision beta,v,root,lnln,nqt,xa,xb
c
      if (v.gt.1.0d3) then
         nqt=2.0d0*((1.0d0-beta)*log(v)-beta*log(beta))
      else
         root=sqrt(v**2+2.0d0*(1.0d0-2.0d0*beta)*v+1.0d0)
c
         lnln=root-(1.0d0-2.0d0*beta+v)
         lnln=lnln/(root+(1.0d0-2.0d0*beta)*v+1.0d0)
         lnln=log(lnln)
         xa=1.0d0-2.0d0*beta
         xb=(xa+v)/root
c
         nqt=xa*log(1.0d0/beta)
         nqt=nqt+xb*lnln
         nqt=nqt+(xb-xa)*log(v)
      endif
c
      fqt=nqt
c      write(6,*)'nqt=',nqt,' v = ',v,' beta = ',beta
      return
      end
c---------------------------------------------------------------
      double precision function fg(z,v)
      implicit DOUBLE PRECISION (a-g,o-z)
      DOUBLE PRECISION ng,lnln

c
c      double precision z,v,root,lnln,ng,xa,xb
c
      if ((1.0d0-z)*v**2/(v-1.0d0)**2.lt.1.0d-3
     +                          .and.v.gt.1.0d0) then
         ng=(3.0d0*v-5.0d0)*v/(v-1.0d0)**2
         ng=ng+2.0d0*(1.0d0+1.0d0/(v-1.0d0)**3)*log(v)
         ng=ng*(1.0d0-z)**2
      elseif (v.gt.1.0d3) then
         ng=2.0d0*(z-1.0d0-z**2*log(z)-(1.0d0-z)**2*log(v))
      else
         root=sqrt(v**2+2.0d0*(1.0d0-2.0d0*z)*v+1.0d0)
c
         lnln=root-(1.0d0-2.0d0*z+v)
         lnln=lnln/(root+(1.0d0-2.0d0*z)*v+1.0d0)
         lnln=log(lnln)
         xa=1.0d0+v-2.0d0*z*(1.0d0-z)
         xb=root-2.0d0*z*(1.0d0-z)/root
c
         ng=xa*log(1.0d0/z)
         ng=ng+(2.0d0*v-xa+xb)*log(v)
         ng=ng+xb*lnln
      endif
      ng=ng/(1.0d0-z)
c       write(6,*)'ng= ',ng,' v = ',v,' beta = ',beta
c
      fg=ng
      return
      end
c---------------------------------------------------------------------
      subroutine xntotal(xnqt,xng,xnql,q2,x_pom,beta,v,z,t)
c
c
      implicit DOUBLE PRECISION (a-g,o-z)
      DOUBLE PRECISION lim
c      double precision xnqt,fqt,xng,fg,xnql,fql,lim,
c     +                 q2,yb,beta,v,z,t,x_pom,jacob,pom,xpom,alp
c
c
      pi=2.0d0*asin(1.0d0)
      q02=1.0d0
      alp=0.25d0
      xng = 0.d0
      xnqt = 0.d0
      xnql = 0.d0
c
c
      xpom=pom(x_pom,-t,v)
      xpom = xpom/x_pom
      if (v.lt.q02/q2) then
         xng=0.0d0
      else
         xng=fg(z,v)/v
         xng=xng**2/q02*((1.0d0-beta/z)**2+(beta/z)**2)
         xng=xng*9.0d0/4.0d0*alp/8.0d0/pi*log(v*q2/q02)
c
         xng = xng/12.d0/z**2 * xpom * beta
c          write(6,*) ' xntotal xng,xpom',xng,xpom

      endif
c
      xnql=fql(beta,v)
      xnql=xnql**2/q2/v**3*4.0d0/3.0d0
c
      xnql = xnql/6.d0 * xpom * beta**3
c
c      write(6,*) ' xntotal xnql,xpom',xnql,xpom
      xnqt=fqt(beta,v)/v
      xnqt=xnqt**2/q02/6.0d0/(1.0d0-beta)
      xnqt=xnqt/sqrt(1.0d0-4.0d0*q02/q2*beta/v)
c add in new factor from ANL-HEP-PR-97-03
cccccccc      xnqt=xnqt*(1.0d0-2.0d0*q02/q2*beta/v)
      xnqt = xnqt/12.d0*beta *xpom
c      write(6,*) ' xntotal xnqt,xpom',xnqt,xpom
c
      lim=1.0d0-4.0d0*q02/q2*beta/v
      if (lim.lt.0.01d0) then
         xnql=xnql*20.0d0
         xnqt=xnqt*20.0d0
      else
         xnql=xnql/dsqrt(lim)
         xnqt=xnqt/dsqrt(lim)
      endif
c
      return
      end
