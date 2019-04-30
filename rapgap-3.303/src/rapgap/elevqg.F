      SUBROUTINE ELEVQG(WT1)
C
C     Matrix elements for 
C     qq --> W/Z g
C     gq --> W/Z q
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
      Integer I,J,NFG,NFLG,NFLP,NFLT,NFLTG,NFLTP,NFLTPC,KPF
      Double Precision draprn
      DOUBLE PRECISION alphasf,gmass
      PARAMETER (alphasf=0.23D0)
      PARAMETER (gmass=.8D0)
      wt1 = 0.d0
      ALPHA_S=ALPHAS(DSQRT(Q2Q))
      if(xpr.GT.1.d0) write(6,*) ' ELEVQG xp2>1',xpr,ipro
      IF(IDIR.EQ.1) THEN
         XPRO = XPR
      ELSE
         XPRO = XPR/XFGKI
      ENDIF
c	write(6,*) ' start ELEVQG',xpro,q2q,t2gki
      if(xpro.ge.0.999) return
      CALL rgstfu(KINT(2,2),XPRO,SNGL(Q2Q),XPQ)
      if(iabs(ke).ne.2212) then
c call virtual photon structure function
         XGAM = XEL/YY
         if(xgam.gt.1.0) then
            write(6,*) ' ELEVQG xpr,x_gam,yy',xpr,xgam,yy,q2q,q2
         endif
         if(xgam.gt.0.999) return
         CALL RYSTGA(XGAM,SNGL(Q2Q),SNGL(Q2),XGQ)
         if(xgq(0).ne.xgq(0)) then
            write(6,*) ' ELEVQG x_gam,Q2Q,Q2,xgq(0)',xgam,q2q,q2,xgq(0)
         endif
      elseif(iabs(ke).eq.2212) then
         if(xel.ge.0.999) return
	   IF(KINT(1,2).NE.990) then
            CALL rgstfu(KINT(1,2),XEL,SNGL(Q2Q),XGQ)
	      ELSE
	      XGQ(0)=1.
	      ENDIF
      else
         write(6,*) ' ELEVQG: KE= ',KE,' not implemented '
         STOP
      endif
c	write(6,*) ' after strfunc'
c	write(6,*) ' ELEVQG KINT ',KINT(1,2),KINT(2,2),xpq(0),xgq(0)
c	call dulist(1)
c      write(6,*) ' ELEVQG ',xgq
c      write(6,*) ' ELEVQG ',NFLAV,IHFLA
      XPQ2 = 0.
      XGQ1 = 0.
      DO 10 I=1,NFLAV
         XPQ2 = XPQ2+XPQ(I) + XPQ(-I)
         XGQ1 = XGQ1+XGQ(I) + XGQ(-I)
   10 CONTINUE
      IF(IHFLA.GE.4) THEN
         XPQ2_HF = XPQ(IHFLA) + XPQ(-IHFLA)
         XGQ1_HF = XGQ(IHFLA) + XGQ(-IHFLA)
      ENDIF
      SH = 2.D0 * AM(1)**2 + 2.D0 * DOT1(NF1,NF2)
      TH = AM(1)**2 - 2.D0 * DOT1(NIA1,NF1) + dot1(nia1,nia1)
      UH = AM(2)**2 - 2.D0 * DOT1(NIA1,NF2) + dot1(nia1,nia1)
c     SUMA qq -> Z/W g
      SUMA = 0 
      SUMA = SUMA * XGQ1*XPQ2
c     SUMB qg -> Z/W q
      SUMB = 1
      SUMB = SUMB * (XGQ(0)*XPQ2 + XGQ1*XPQ(0))
      

c now the total
      SUTO = SUMA + SUMB 
      SUM2 = SUTO*16.D0*PI*PI*ALPHA_S*ALPHA_S
      WT1=SUM2
      RNTEST = draprn()
      IF((SUMA)/SUTO.LT.RNTEST) THEN
         IRESPRO = 2
      ELSE
         IRESPRO = 1
      ENDIF
      ICOLORA = 0
c      write(6,*) ' ELEVQG ',IRESPRO
      IF(IRESPRO.EQ.1) THEN
c qq -> Z/W g
c select flavour of incoming q's
         Qtot = 0.
         Do I=1,NFLAV
           Qtot = Qtot + XGQ(I)*XPQ(-I) +  XGQ(-I)*XPQ(I)
         End do
         NFLT = -NFLAV-1
         QSUM = -QTOT*draprn()
   50    NFLT = NFLT + 1
         IF(NFLT.EQ.0) GOTO 50
            QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(-NFLT))
            IF(QSUM.LE.0.0D0) GOTO 50
         IF(NFLT.GT.NFLAV) write(6,*) ' ELEVQG NFL > NFLAV ',NFLT,
     +      NFLAV
         K(NIA1,2) = -NFLT    
         K(NIR1,2) = NFLT
         K(NIA2,2) = NFLT
c select flavour of outgoing q q_bar
         K(NF2,2) = 21
      Elseif(IRESPRO.EQ.2) THEN
c     SUMB qg -> Z/W q
c     quark from hadron 1
         NFLT = -NFLAV-1
c         NFLT =0
         QSUM = -DBLE(XGQ1)*draprn()
  51     NFLT = NFLT + 1
         IF(NFLT.EQ.0) GOTO 51
c               write(6,*) ' test ',QSUM,XPQC1,NFLT,XPQ2(NFLT)
         QSUM = QSUM + DBLE(XGQ(NFLT))
         IF(QSUM.LE.0.0D0) GOTO 51
         IF(NFLT.GT.NFLAV) write(6,*) ' meoffsh NFL > NFLAV ',
     +      NFLT, NFLAV
         K(NIA2,2) = 21
         K(NIA1,2) = NFLT
         K(NF2,2) =  NFLT
      ENDIF
      if(NFLT.EQ.0) then 
      write(6,*) ' NFLT ',NFLT
      endif
      SCAL1 = SNGL(Q2Q)    	
      if(K(NIA1,2).EQ.21) then
         XPD1 = XGQ(0)
	   ELSE
         XPD1 = XGQ(K(NIA1,2))
	   ENDIF	
c	write(6,*) ' end of ELEVQG',wt1
 
      IF(WT1.LE.0.D0.AND.IDEBUG.EQ.1) THEN
         write(6,*) ' ELEVQG wt1 = ',wt1,' Q2 = ',Q2,' Q2Q = ',Q2Q
         write(6,*) ' ELEVQG xgq = ',(xgq(i),i=-6,6)
         write(6,*) ' ELEVQG xpq = ',(xpq(i),i=-6,6)
      ENDIF
      IF(WT1.LT.0.D0) THEN
         write(6,*) ' ELEVQG wt1 = ',wt1,' Q2 = ',Q2,' Q2Q = ',Q2Q
         WT1 = 0.D0
      ENDIF
      RETURN
      END
