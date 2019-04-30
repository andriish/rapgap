      SUBROUTINE F2MCD(BETA,X_POM,Q2,T2,XDX,FLQ)
* calculate   in the approach of M. McDermott
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
#include "rgparas.inc"
#include "rgparam.inc"
#include "rgluco.inc"
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
      COMMON /SEL/IPHI
      DOUBLE PRECISION XDX(-4:4)
      DOUBLE PRECISION C1,Cg
      COMMON/BUCHMUE/C1,Cg
#include "rgludat2.inc"
      EXTERNAL FMD
      IPHI = 0
c      C1 = 1.d0
c      Cg = 1.d0
      q2t = q2
      t2t = t2
      bet = beta
      xpom = x_pom
      ymcd = xpom*bet
      umcd = xpom/ymcd - 1.d0
c      write(6,*) ' xpom',x_pom,beta,q2,t2
c include Donnacchi/Landshoff dipole form factor for t
      dl =(4.d0-2.8d0*t2t)/(4.d0-t2t)*
     +   1.d0/(1.d0-t2t/0.7d0)**2
      DO 10 I=-4,4
   10 XDX(I)=0.d0
c 1/xpom removed due to jacobian dy/dx_g=x_pom
      XDX(0) = C1/(cg + umcd)/xpom/bet*dl
      XDX(-4)=XDX(4)
      XDX(-3)=XDX(3)
      XDX(-2)=XDX(2)
      XDX(-1)=XDX(1)
      RETURN
      END

