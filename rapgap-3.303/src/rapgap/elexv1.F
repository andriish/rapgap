      SUBROUTINE ELEXV1(WT1)
C
C     Matrix elements for qq -> W/Z
C
C
      IMPLICIT NONE
#include "rglujets.inc"
#include "rgpara1.inc"
#include "rgpara.inc"
#include "rgluco.inc"
#include "rgpart.inc"
#include "rgrapgki.inc"
#include "rgludat2.inc"
#include "rgstru.inc"
#include "rgdisdif.inc"
#include "rgcolcon.inc"
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG

      Integer IRFUL
      REAL XPQ(-6:6),XGQ(-6:6)
      REAL XGQ1,XPQ2,XGAM,XPRO,XGQ1_HF,XPQ2_HF
      REAL SNGL
c IRFUL=1 using full ME's
c IRFUL=0 using ME's from Bengtsson
      PARAMETER (IRFUL=0)
ccc      PARAMETER (IRFUL=1)
      DOUBLE PRECISION ZETA3
C---ZETA3=RIEMANN ZETA FUNCTION(3)
      PARAMETER (ZETA3=1.202056903159594D0)
C---paramtric fit to the function used in the Bartels/Forshaw/et al.
calculation
      Double Precision wt1,ALPHA_S,ALPHAS,SH,UH,TH,VTH,VUH
      Double Precision SUMA,SUMB,SUMC,SUMC1,SUMC1T,SUMD,SUMD1,SUMD2
      Double Precision SUMF,SUME,SUME1,SUME2,SUMTE
      Double Precision SUMTF,SUMF1,SUMF2,SUMG,SUMG1,SUTO,SUM2
      Double Precision A1,A2,A3,ASUM,FLSUM,FLRN
      Double Precision QTOT,QSUM,CF,CA,yh,omeg0,delta
      Double Precision RNTEST
      Double Precision PDF
      Integer I,J,NFLT
      Double Precision draprn
      DOUBLE PRECISION alphasf,gmass
      Double Precision gDY,GF,Vq2,Aq2,M2,Qf,Tf
      PARAMETER (alphasf=0.23D0)
      PARAMETER (gmass=.8D0)
      wt1 = 0.d0
      M2 = AM(1)**2
      ALPHA_S=ALPHAS(DSQRT(Q2Q))
      if(xpr.GT.1.d0) write(6,*) ' ELEXV1 xp2>1',xpr,ipro
      IF(IDIR.EQ.1) THEN
         XPRO = XPR
      ELSE
         XPRO = XPR/XFGKI
      ENDIF
c	write(6,*) ' start ELEXV1',xpro,q2q,t2gki
      if(xpro.ge.0.999) return
      CALL rgstfu(KINT(2,2),XPRO,SNGL(Q2Q),XPQ)
      if(iabs(ke).ne.2212) then
c call virtual photon structure function
         XGAM = XEL/YY
         if(xgam.gt.1.0) then
            write(6,*) ' ELEXV1 xpr,x_gam,yy',xpr,xgam,yy,q2q,q2
         endif
         if(xgam.gt.0.999) return
         CALL RYSTGA(XGAM,SNGL(Q2Q),SNGL(Q2),XGQ)
         if(xgq(0).ne.xgq(0)) then
            write(6,*) ' ELEXV1 x_gam,Q2Q,Q2,xgq(0)',xgam,q2q,q2,xgq(0)
         endif
      elseif(iabs(ke).eq.2212) then
         if(xel.ge.0.999) return
	   IF(KINT(1,2).NE.990) then
            CALL rgstfu(KINT(1,2),XEL,SNGL(Q2Q),XGQ)
	      ELSE
	      XGQ(0)=1.
	      ENDIF
      else
         write(6,*) ' ELEXV1: KE= ',KE,' not implemented '
         STOP
      endif
c	write(6,*) ' after strfunc'
c	write(6,*) ' ELEXV1 KINT ',KINT(1,2),KINT(2,2),xpq(0),xgq(0)
c	call dulist(1)
c      write(6,*) ' ELEXV1 ',xgq
c      write(6,*) ' ELEXV1 ',NFLAV,IHFLA
      XPQ2 = 0.
      XGQ1 = 0.
      DO 10 I=1,NFLAV
         XPQ2 = XPQ2+XPQ(I) + XPQ(-I)
         XGQ1 = XGQ1+XGQ(I) + XGQ(-I)
   10 CONTINUE


      SUMA = 1.D0
c now the total
c calculate it for gamma production 
c        16\pi^2*alpha_s*\alpha_em *ME^2         
c calculate it for W production (from Ellis Stirling Webber, eq 9.53 +corr)
c        4\pi*alpha_s*\sqrt{2}*GF*Mw^2*|Vqq|^2 *ME^2         
c calculate it for Z production (from Ellis Stirling Webber, eq 9.53 +corr)
c        4\pi*alpha_s*\sqrt{2}*GF*Mz^2*(Vq^2 + Aq^2) *ME^2 
c        gDY =  coupling to vector boson
c        dDY for W : \pi * g_weak *|Vqq|^2 
c        dDY for Z0: \pi * g_weak *(Vq^2 + Aq^2) 
c        g_weak = sqrt(2)*GF *Mz^2*(Vq^2 + Aq^2)
	   GF = 1.166E-5
c  vector and axial couplings for Z         
c                (u,c,t, 2/3)                  (d,s,b, -1/3)         
c         Vq2 = 2*(0.5 - 2.*2./3.*sin2w)**2 + 3*(-0.5 + 2.*1./3.*sin2w)**2
c         Aq2 = 2*0.5*0.5 + 3*0.5*0.5
         gDY = 4.*pi *sqrt(2.)* M2 * GF  
c     SUMA qq -> Z/W g
      SUMA = 1./3./4.
      PDF = 0
      DO I=1,NFLAV
         If(I.eq.1.or.I.eq.3.or.I.eq.5) then
            Qf = -1./3.
            Tf = 1./2.
         else 
            Qf = 2./3.
            Tf = -1./2.
         Endif
c         write(6,*) i,Qf,Tf
         Vq2 = (Tf - 2.*Qf*sin2w)**2
         Aq2 = Tf**2
         PDF = PDF + (Vq2 + Aq2)*(XPQ(I)* XGQ(-I) + XPQ(-I)* XGQ(I))
      End do
c      write(6,*) ' suma, gDY,PDF', suma, gDY,PDF
      SUMA = SUMA * PDF

      SUM2 = SUMA*gDY

c now the total
      WT1=SUM2
      Qtot = 0
         Do I=1,NFLAV
           Qtot = Qtot + XGQ(I)*XPQ(-I) +  XGQ(-I)*XPQ(I)
         End do
         If(Qtot.eq.0) then
            WT1 = 0.D0
            Return
         Endif

         NFLT = -NFLAV-1
         QSUM = -QTOT*draprn()
   50    NFLT = NFLT + 1
         IF(NFLT.EQ.0) GOTO 50
          IF(NFLT.GT.NFLAV) write(6,*) ' ELEVQG NFL > NFLAV ',NFLT,
     +      NFLAV,Qtot
           QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(-NFLT))
c            write(6,*) nflt,Qsum,DBLE(XPQ(NFLT)*XGQ(-NFLT))
            IF(QSUM.LE.0.0D0) GOTO 50
         K(NIA1,2) = -NFLT    
         K(NIR1,2) = NFLT
         K(NIA2,2) = NFLT
c select flavour of outgoing parton
         K(NF2,2) = 21

c	write(6,*) ' end of ELEXV1',wt1

      IF(WT1.LE.0.D0.AND.IDEBUG.EQ.1) THEN
         write(6,*) ' ELEXV1 wt1 = ',wt1,' Q2 = ',Q2,' Q2Q = ',Q2Q
         write(6,*) ' ELEXV1 xgq = ',(xgq(i),i=-6,6)
         write(6,*) ' ELEXV1 xpq = ',(xpq(i),i=-6,6)
      ENDIF
      IF(WT1.LT.0.D0) THEN
         write(6,*) ' ELEXV1 wt1 = ',wt1,' Q2 = ',Q2,' Q2Q = ',Q2Q
         WT1 = 0.D0
      ENDIF
      RETURN
      END
