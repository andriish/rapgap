      SUBROUTINE ELEH0(WT1)
C
C     Matrix elements for resolved photon processes
C     gg --> qq_bar
C     g + g --> g + g
C     g + q --> g + q
C     qq_bar --> gg
C     q + q_bar --> q + q_bar
C     qq --> qq
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
      if(xpr.GT.1.d0) write(6,*) ' ELEH0 xp2>1',xpr,ipro
      IF(IDIR.EQ.1) THEN
         XPRO = XPR
      ELSE
         XPRO = XPR/XFGKI
      ENDIF
c	write(6,*) ' start ELEH0',xpro,q2q,t2gki
      if(xpro.ge.0.999) return
      CALL rgstfu(KINT(2,2),XPRO,SNGL(Q2Q),XPQ)
      if(iabs(ke).ne.2212) then
c call virtual photon structure function
         XGAM = XEL/YY
         if(xgam.gt.1.0) then
            write(6,*) ' ELEH0 xpr,x_gam,yy',xpr,xgam,yy,q2q,q2
         endif
         if(xgam.gt.0.999) return
         CALL RYSTGA(XGAM,SNGL(Q2Q),SNGL(Q2),XGQ)
         if(xgq(0).ne.xgq(0)) then
            write(6,*) ' ELEH0 x_gam,Q2Q,Q2,xgq(0)',xgam,q2q,q2,xgq(0)
         endif
      elseif(iabs(ke).eq.2212) then
         if(xel.ge.0.999) return
	   IF(KINT(1,2).NE.990) then
            CALL rgstfu(KINT(1,2),XEL,SNGL(Q2Q),XGQ)
	      ELSE
	      XGQ(0)=1.
	      ENDIF
      else
         write(6,*) ' ELEH0: KE= ',KE,' not implemented '
         STOP
      endif
c	write(6,*) ' after strfunc'
c	write(6,*) ' ELEH0 KINT ',KINT(1,2),KINT(2,2),xpq(0),xgq(0)
c	call dulist(1)
c      write(6,*) ' ELEH0 ',xgq
c      write(6,*) ' ELEH0 ',NFLAV,IHFLA
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
      VTH = - 2.D0 * DOT1(NIA1,NF1)
      VUH = - 2.D0 * DOT1(NIA1,NF2)
      SUMA = 0.D0
      SUMB = 0.D0
      SUMC = 0.D0
      SUMD = 0.D0
      SUME = 0.D0
      SUMF = 0.D0
      SUMG = 0.D0
      IF(IRPA.NE.0) THEN
         IF(IHFLA.LT.4) THEN
C.....MATRIX ELEMENT for gg --> qq_barp(466)
            IF(IRFUL.EQ.1) THEN
               SUMA = 2.d0*th*uh/th**2/12.d0 + 2.d0*th*uh/uh**2/12.d0 -
     +         3.d0/16.d0*4.d0*(1.d0 -uh*th/sh**2) + 3.d0/32.d0*4.d0
            ELSE
c this was from bengtsson
               SUMA = uh/th -2.d0*uh**2/sh**2 + th/uh - 2.d0*th**2/sh**
     +         2
               SUMA = SUMA/6.d0
            ENDIF
c include sum over NFG flavors and parton densities
            NFG = 3
            IF(SH.GT.4.D0*DBLE(PMAS(4,1)**2)) NFG = 4
            IF(SH.GT.4.D0*DBLE(PMAS(5,1)**2)) NFG = 5

            SUMA = SUMA * DFLOAT(NFG) * DBLE(XGQ(0)) * DBLE(XPQ(0))
         ELSEIF(IHFLA.GE.4) THEN
c heavy flavour gg --> QQ_bar
            SUMA = (32.D0/vth/vuh - 8d0*9D0/sh**2)*
     +   (vth**2/4d0 + vuh**2/4d0 + AM(1)**2*sh - AM(1)**2*sh/vth/vuh)
            SUMA = SUMA/16d0/3d0
            SUMA = SUMA  * DBLE(XGQ(0)) * DBLE(XPQ(0))
         ENDIF
      ENDIF
      IF(IRPB.NE.0) THEN
         IF(IRFUL.EQ.1) THEN
C.....MATRIX ELEMENT for gg --> gg
            SUMB = 9.d0/8.d0*(17.d0/2.d0 - 4.d0*uh*sh/th**2 + 17.d0/
     +      2.d0 - 4.d0*sh*th/uh**2 + 17.d0/2.d0 - 4.d0*uh*th/sh**2 +
     +      27.d0) + 9.d0/16.d0*(15.d0 - sh**2/th/uh + 15.d0 - uh**2/
     +      th/sh - 15.d0 + th**2/sh/uh) - 3.d0*9.d0/8.d0*81.d0/4.d0
         ELSE
c this is from bengtsson
            SUMB = 9.d0/4.d0*( sh**2/th**2+2.d0*sh/th+3.d0+2.d0*th/sh+
     +      th**2/sh**2 +uh**2/sh**2+2.d0*uh/sh+3.d0+2.d0*sh/uh+sh**2/
     +      uh**2 +th**2/uh**2+2.d0*th/uh+3.d0+2.d0*uh/th+uh**2/th**2 )
     +      * 0.5d0
ccc         write(6,*) sumb,xgq(0),xpq(0),uh,th,sh
         ENDIF
         SUMB = SUMB * DBLE(XGQ(0)) * DBLE(XPQ(0))

      ENDIF
      IF(IRPC.NE.0) THEN
C.....MATRIX ELEMENT for qg --> qg
         IF(IRFUL.EQ.1) THEN
            SUMC1 = 2.d0*(1.d0-uh*sh/th**2)-4.d0/9.d0*(sh/uh+uh/sh)-
     +      1.d0
c  this was from bengtsson
         ELSE
            SUMC1 = 4.d0/9.d0*(2.d0*uh**2/th**2-uh/sh+2.d0*sh**2/th**2-
     +      sh/ uh)
         ENDIF
         sumc1t = SUMC1 *(DBLE(XGQ1) * DBLE(XPQ(0)))
c      write(6,*) xgq1,xgg1,xpq2,xpg2
         SUMC = SUMC1 *(DBLE(XGQ1)*DBLE(XPQ(0))+DBLE(XGQ(0))*DBLE(XPQ2)
     +   )
         IF(IHFLA.GE.4) THEN
            SUMC=0.D0
            IF(sh.ge.4.*PYMASS(IHFLA)) THEN
               SUMC = SUMC1* (DBLE(XGQ1_HF)*DBLE(XPQ(0))+DBLE(XGQ(0))*
     +         DBLE(XPQ2_HF))
            ENDIF
         ENDIF
      ENDIF
      IF(IRPD.NE.0) THEN
C.....MATRIX ELEMENT for qq_bar --> gg
         IF(IRFUL.EQ.1) THEN
            SUMD1 = 16.d0/27.d0*2.d0*uh/th + 16.d0/27.d0*2d0*th/uh
     +      -4.d0/3.d0*4.d0*(1.d0-uh*th/sh**2) +2.d0/3.d0*4.d0
            SUMD2 = 0.D0
         ELSE
            SUMD1 =32.d0/27.d0*(uh/th-2.d0*uh**2/sh**2)
            SUMD2 =32.d0/27.d0*(th/uh-2.d0*th**2/sh**2)
         ENDIF
c
         PDF = 0
         Do I=1,NFLAV
         PDF = PDF + XGQ(I)*XPQ(-I)
         ENDDO
         SUMD = (SUMD1+SUMD2) * PDF
      ENDIF
      IF(IRPE.NE.0) THEN
C.....MATRIX ELEMENT for qq_bar --> qq_bar
c here take only s-channel diagram (as done in PYTHIA)

         IF(IHFLA.LT.4) THEN

            SUME1 = 4.d0/9.d0*(sh**2 + uh**2)/th**2
            SUME2 = 4.d0/9.d0*(th**2 + uh**2)/sh**2
c include sum over NFG flavors and parton densities
            NFG = 3
            IF(SH.GT.4.D0*DBLE(PMAS(4,1)**2)) NFG = 4
            IF(SH.GT.4.D0*DBLE(PMAS(5,1)**2)) NFG = 5
            SUME2 = SUME2*NFG
         ELSEIF(IHFLA.GE.4) THEN
c heavy flavour qq --> QQ_bar
            IF(sh.ge.4.*PYMASS(IHFLA)) THEN
            SUME1 = 4.d0/9.d0*(sh**2 + uh**2)/th**2
            SUME2=4.d0/9.d0*
     +      ((vth**2 + vuh**2)/sh**2 + 4d0*AM(1)**2/sh**2)
            endif
         ENDIF
c
c first take t channel: q_i + q_bar_i --> q_i + q_bar_i
         IF(IHFLA.LT.4) THEN
         PDF = 0.
         Do I=1,NFLAV
           Do J=1,NFLAV
             PDF = PDF + XGQ(I)*XPQ(-J) + XGQ(-I)*XPQ(J) 
           end do 
         end do
         SUME1 = SUME1 * PDF
         ELSEIF(IHFLA.EQ.4) THEN
         SUME1 = SUME1 * (
     +                    DBLE(XGQ(4))*DBLE(XPQ(-4)) +
     +			  DBLE(XGQ(-4))*DBLE(XPQ(4)))
         ELSEIF(IHFLA.EQ.5) THEN
         SUME1 = SUME1 * (
     +                    DBLE(XGQ(5))*DBLE(XPQ(-5)) +
     +			  DBLE(XGQ(-5))*DBLE(XPQ(5)))

         ENDIF
c now s channel q_i + q_bar_i --> q_k + q_bar_k
         PDF = 0.
         Do I=1,NFLAV
            PDF = PDF + XGQ(I)*XPQ(-I) + XGQ(-I)*XPQ(I)
         end do
         SUME2 = SUME2 * PDF
         SUME = SUME1 + SUME2
      ENDIF
      IF(IRPF.NE.0) THEN
C.....MATRIX ELEMENT for qq --> qq, q_bar q_bar --> q_bar q_bar


         SUMF1 = 4.d0/9.d0*(sh**2 + uh**2)/th**2
         SUMF2 = 4.d0/9.d0*(sh**2 + th**2)/uh**2
c

         IF(IHFLA.LT.4) THEN
         PDF = 0.
         Do I=1,NFLAV
           Do J=1,NFLAV
             PDF = PDF + XGQ(I)*XPQ(J) + XGQ(-I)*XPQ(-J) 
           end do 
         end do
         SUMF1 = SUMF1 * PDF
         ELSEIF(IHFLA.EQ.4) THEN
         SUMF1 = SUMF1 * (DBLE(XGQ(1))*DBLE(XPQ(4)) + 
     +                    DBLE(XGQ(2))*DBLE(XPQ(4)) +
     +                    DBLE(XGQ(3))*DBLE(XPQ(4)) +
     +                    DBLE(XGQ(4))*DBLE(XPQ(4)+XPQ(3)+XPQ(2)+
     +                    XPQ(1)) +
     +                    DBLE(XGQ(-1))*DBLE(XPQ(-4)) + 
     +                    DBLE(XGQ(-2))*DBLE(XPQ(-4)) + 
     +		        DBLE(XGQ(-3))*DBLE(XPQ(-4)) + 
     +			  DBLE(XGQ(-4))*DBLE(XPQ(-4)+XPQ(-3)+XPQ(-2)+
     +                    XPQ(1))+
     +                    DBLE(XGQ(4))*DBLE(XPQ(-4)) +
     +			  DBLE(XGQ(-4))*DBLE(XPQ(4)))
         ELSEIF(IHFLA.EQ.5) THEN
         SUMF1 = SUMF1 * (DBLE(XGQ(1))*DBLE(XPQ(5)) + 
     +                    DBLE(XGQ(2))*DBLE(XPQ(5)) +
     +                    DBLE(XGQ(3))*DBLE(XPQ(5)) +
     +                    DBLE(XGQ(4))*DBLE(XPQ(5)) +
     +                    DBLE(XGQ(5))*DBLE(XPQ(5)+XPQ(4)+XPQ(3)+
     +                    XPQ(2)+XPQ(1)) +
     +                    DBLE(XGQ(-1))*DBLE(XPQ(-5)) + 
     +                    DBLE(XGQ(-2))*DBLE(XPQ(-5)) + 
     +		        DBLE(XGQ(-3))*DBLE(XPQ(-5)) + 
     +		        DBLE(XGQ(-4))*DBLE(XPQ(-5)) + 
     +			  DBLE(XGQ(-5))*DBLE(XPQ(-5)+XPQ(-4)+XPQ(-3)+
     +                    XPQ(-2)+XPQ(-1)) +
     +                    DBLE(XGQ(5))*DBLE(XPQ(-5)) +
     +			  DBLE(XGQ(-5))*DBLE(XPQ(5)))
     

         ENDIF
         IF(IHFLA.LT.4) THEN
         PDF = 0.
         Do I=1,NFLAV
            PDF = PDF + XGQ(I)*XPQ(I) + XGQ(-I)*XPQ(-I)
         end do
         SUMF2 = SUMF2 *PDF
         ELSEIF(IHFLA.EQ.4) THEN
         SUMF2 = SUMF2 * (DBLE(XGQ(1))*DBLE(XPQ(4)) + 
     +                    DBLE(XGQ(2))*DBLE(XPQ(4)) +
     +                    DBLE(XGQ(3))*DBLE(XPQ(4)) +
     +                    DBLE(XGQ(4))*DBLE(XPQ(4)+XPQ(3)+XPQ(2)+
     +                    XPQ(1)) +
     +                    DBLE(XGQ(-1))*DBLE(XPQ(-4)) + 
     +                    DBLE(XGQ(-2))*DBLE(XPQ(-4)) + 
     +		        DBLE(XGQ(-3))*DBLE(XPQ(-4)) + 
     +			  DBLE(XGQ(-4))*DBLE(XPQ(-4)+XPQ(-3)+XPQ(-2)+
     +                    XPQ(1)))
         ELSEIF(IHFLA.EQ.5) THEN
         SUMF2 = SUMF2 * (DBLE(XGQ(1))*DBLE(XPQ(5)) + 
     +                    DBLE(XGQ(2))*DBLE(XPQ(5)) +
     +                    DBLE(XGQ(3))*DBLE(XPQ(5)) +
     +                    DBLE(XGQ(4))*DBLE(XPQ(5)) +
     +                    DBLE(XGQ(5))*DBLE(XPQ(5)+XPQ(4)+XPQ(3)+
     +                    XPQ(2)+XPQ(1)) +
     +                    DBLE(XGQ(-1))*DBLE(XPQ(-5)) + 
     +                    DBLE(XGQ(-2))*DBLE(XPQ(-5)) + 
     +		        DBLE(XGQ(-3))*DBLE(XPQ(-5)) + 
     +		        DBLE(XGQ(-4))*DBLE(XPQ(-5)) + 
     +			  DBLE(XGQ(-5))*DBLE(XPQ(-5)+XPQ(-4)+XPQ(-3)+
     +                    XPQ(-2)+XPQ(-1)))
         ENDIF
         SUMF = SUMF1 + SUMF2
      ENDIF
      IF(IRPG.NE.0) THEN
c cross sections from B. Cox J. Forshaw hep-ph/9805206
         SUMG1 = 0.D0
         if(IRPG.eq.1) then
c eq. 3 Mueller tang
            CF = 4.D0/3.D0
            CA = 3.D0
            yh=LOG(SH/(-TH))+1.D0
            omeg0 = CA*alphasf/pi*4*dlog(2.d0)
            SUMG1=16.D0*PI*2*PI**3*(CF*alphasf)**4
     +    *(SH/TH)**2 *dexp(2.d0*omeg0*yh)/(7.D0*alphasf*CA*zeta3*yh)**3
         elseif(IRPG.eq.2) then
c eq 4. massive gluon
            delta=sqrt(1-((4*gmass**2)/TH))
            SUMG1=(SH/TH)**2/delta**2*(log((delta+1)/(delta-1)))**2*
     +      PI**2*(alphasf)**4
         endif
         SUMG=SUMG1*DBLE(XPQ2)*DBLE(XGQ1)
      ENDIF
c now the total
      SUTO = SUMA + SUMB + SUMC + SUMD + SUME + SUMF + SUMG
c      SUTO = SUMA + SUMB + SUMC + SUMD + SUME + SUMF
      SUM2 = SUTO*16.D0*PI*PI*ALPHA_S*ALPHA_S
      WT1=SUM2
      RNTEST = draprn()
      IF((SUMA+SUMB+SUMC+SUMD+SUME+SUMF)/SUTO.LT.RNTEST) THEN
         IRESPRO = 7
      ELSEIF((SUMA+SUMB+SUMC+SUMD+SUME)/SUTO.LT.RNTEST) THEN
         IRESPRO = 6
      ELSEIF((SUMA+SUMB+SUMC+SUMD)/SUTO.LT.RNTEST) THEN
         IRESPRO = 5
      ELSEIF((SUMA+SUMB+SUMC)/SUTO.LT.RNTEST) THEN
         IRESPRO = 4
      ELSEIF((SUMA+SUMB)/SUTO.LT.RNTEST) THEN
         IRESPRO = 3
      ELSEIF(SUMA/SUTO.LT.RNTEST) THEN
         IRESPRO = 2
      ELSE
         IRESPRO = 1
      ENDIF
      ICOLORA = 0
c      write(6,*) ' ELEH0 ',IRESPRO
      IF(IRESPRO.EQ.1) THEN
c get color configuartion for gg --> qq_bar
         A1 = uh/th - 2.d0*uh**2/sh**2
         A2 = th/uh - 2.d0*th**2/uh**2
         ASUM = A1 + A2
         IF(A1/ASUM.lt.draprn()) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
         IF(IHFLA.LT.4) THEN
c select flavour of outgoing q q_bar
            FLRN = draprn()
            FLSUM = 1.D0/3.D0
            IF(FLRN.LT.FLSUM) THEN
               KPF = 1
            ELSEIF(FLRN.LT.2.D0*FLSUM) THEN
               KPF = 2
            ELSEIF(FLRN.LE.3.D0*FLSUM) THEN
               KPF = 3
            ENDIF
         ELSE
            KPF = IHFLA
         ENDIF
         IF(draprn().LE.0.5) THEN
            K(NF1,2) = KPF
            K(NF2,2) = - KPF
         ELSE
            K(NF1,2) = - KPF
            K(NF2,2) = KPF
         ENDIF

      ELSEIF(IRESPRO.EQ.2) THEN
c get color configuartion for gg --> gg
         K(NF1,2)=21
         K(NF2,2)=21
         A1 = sh**2/th**2+2.d0*sh/th+3.d0+2.d0*th/sh+th**2/sh**2
         A2 = uh**2/sh**2+2.d0*uh/sh+3.d0+2.d0*sh/uh+sh**2/uh**2
         A3 = th**2/uh**2+2.d0*th/uh+3.d0+2.d0*uh/th+uh**2/th**2
         ASUM = A1+A2+A3
         RNTEST = draprn()
         ICOLORA = 0
         IF((A1+A2)/ASUM.lt.rntest) THEN
            ICOLORA=3
         ELSEIF(A1/ASUM.lt.rntest) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
      ELSEIF(IRESPRO.EQ.3) THEN
c where is the quark ( from gamma or proton)
         IF(SUMC1t/SUMC.GT.draprn()) THEN
            IF(IHFLA.LT.4) THEN
               NFLT = -NFLAV-1
               QSUM = -DBLE(XGQ1)*draprn()
   20          NFLT = NFLT + 1
               IF(NFLT.EQ.0) GOTO 20
c         write(6,*) ' q of gamma ',QSUM,XGQ1,NFLT,XGQ(NFLT)
               QSUM = QSUM + DBLE(XGQ(NFLT))
               IF(QSUM.LE.0.0D0) GOTO 20
               IF(NFLT.GT.NFLAV) write(6,*) ' ELEH0 NFL > NFLAV ',
     +         NFLT, NFLAV
            ELSE
               NFLT = IHFLA
            ENDIF
            IF(draprn().GT.0.5) NFLT=-IHFLA
            K(NIA1,2) = NFLT
            K(NIR1,2) = -NFLT
            K(NIA2,2) = 21
            K(NF1,2) = NFLT
            K(NF2,2) = 21
         ELSE
            IF(IHFLA.LT.4) THEN
               NFLT = -NFLAV-1
               QSUM = -DBLE(XPQ2)*draprn()
   30          NFLT = NFLT + 1
               IF(NFLT.EQ.0) GOTO 30
c         write(6,*) ' q of proton ',QSUM,XPQ2,NFLT,XPQ(NFLT)
               QSUM = QSUM + DBLE(XPQ(NFLT))
               IF(QSUM.LE.0.0D0) GOTO 30
               IF(NFLT.GT.NFLAV) write(6,*) ' ELEH0 NFL > NFLAV ',
     +         NFLT, NFLAV
            ELSE
               NFLT = IHFLA
            ENDIF
            IF(draprn().GT.0.5) NFLT=-IHFLA
            K(NIA1,2) = 21
            K(NIR1,2) = 21
            K(NIA2,2) = NFLT
            K(NF1,2) = 21
            K(NF2,2) = NFLT
         ENDIF
c         call dulist(1)
c get color configuartion for gq --> gq
         A1 = 2.d0*uh**2/th**2 - uh/sh
         A2 = 2.d0*sh**2/th**2 - sh/th
         ASUM = A1 + A2
         ICOLORA = 0
         IF(A1/ASUM.lt.draprn()) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
      ELSEIF(IRESPRO.EQ.4) THEN
c get color configuartion for qq_bar --> gg
         A1 = uh/th - 2.d0*uh**2/sh**2
         A2 = th/uh - 2.d0*th**2/uh**2
         ASUM = A1 + A2
         IF(A1/ASUM.lt.draprn()) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
         IF(IHFLA.LT.4) THEN
            Qtot = 0
            Do I=1,NFLAV
              Qtot = Qtot + XGQ(I)*XPQ(-I)+ XGQ(-I)*XPQ(I)
            End do
            NFLT = -NFLAV-1
            QSUM = -QTOT*draprn()
   40       NFLT = NFLT + 1
            IF(NFLT.EQ.0) GOTO 40
            QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(-NFLT))
            IF(QSUM.LE.0.0D0) GOTO 40
            IF(NFLT.GT.NFLAV) write(6,*) ' ELEH0 NFL > NFLAV ',NFLT,
     +      NFLAV
         ELSE
            NFLT = IHFLA
            IF(draprn().GT.0.5) NFLT= -IHFLA
         ENDIF
c         write(6,*) ' IRPE 2nd ',NFLT
         K(NIA1,2) = -NFLT
         K(NIR1,2) = NFLT
         K(NIA2,2) = NFLT
         K(NF1,2) = 21
         K(NF2,2) = 21
      ELSEIF(IRESPRO.EQ.5) THEN
c qq_bar --> qq_bar
         SUMTE = SUME1+SUME2
         IF(SUME1/SUMTE.lt.draprn()) THEN
            ICOLORA = 2
            Qtot = 0.
            Do I=1,NFLAV
              Qtot = Qtot + XGQ(I)*XPQ(-I) +  XGQ(-I)*XPQ(I)
            End do
            NFLT = -NFLAV-1
            QSUM = -QTOT*draprn()
   50       NFLT = NFLT + 1
            IF(NFLT.EQ.0) GOTO 50
            QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(-NFLT))
            IF(QSUM.LE.0.0D0) GOTO 50
            IF(NFLT.GT.NFLAV) write(6,*) ' ELEH0 NFL > NFLAV ',NFLT,
     +      NFLAV
            K(NIA1,2) = -NFLT
            K(NIR1,2) = NFLT
            K(NIA2,2) = NFLT
c select flavour of outgoing q q_bar
            IF(IHFLA.LT.4) THEN
               FLRN = draprn()
               FLSUM = 1.D0/3.D0
               IF(FLRN.LT.FLSUM) THEN
                  KPF = 1
               ELSEIF(FLRN.LT.2.D0*FLSUM) THEN
                  KPF = 2
               ELSEIF(FLRN.LE.3.D0*FLSUM) THEN
                  KPF = 3
               ENDIF
            ELSE
               KPF = IHFLA
            ENDIF
            IF(K(NIA1,2).LE.0) THEN
               K(NF1,2) = -KPF
               K(NF2,2) =  KPF
            ELSE
               K(NF1,2) = KPF
               K(NF2,2) = -KPF
            ENDIF
         ELSE
            ICOLORA = 1
		IF(IHFLA.LT.4) THEN
               Qtot = 0.
               Do I=1,NFLAV
                 Do J=1,NFLAV
                    Qtot = Qtot + XGQ(I)*XPQ(-J) + XGQ(-I)*XPQ(J)
                 end do
               end do
               NFLTG = -NFLAV-1
               QSUM = -QTOT*draprn()
   60          NFLTG = NFLTG + 1
               IF(NFLTG.EQ.0) GOTO 60
               IF(NFLTG.LE.0) THEN
                  NFLTP = 0
                  NFLTPC = NFLAV
               ELSE
                  NFLTP = -NFLAV-1
                  NFLTPC = 0
               ENDIF
   70          NFLTP = NFLTP + 1
               IF(NFLTP.EQ.0) GOTO 70
               IF(NFLTP.GT.NFLTPC) GOTO 60
               QSUM = QSUM + DBLE(XPQ(NFLTP)*XGQ(NFLTG))
               IF(QSUM.LE.0.D0.AND.NFLTP.LE.NFLTPC) GOTO 70
               IF(QSUM.LE.0.0D0) GOTO 60
		ELSE
		   QTOT =DBLE(
     +	         (XGQ(1)+XGQ(2)+XGQ(2)+XGQ(3)+XGQ(4))*XPQ(-4))+
     +               DBLE(XGQ(4)*(XPQ(-4)+XPQ(-3)+XPQ(-2)+XPQ(-1)))+
     + 	         DBLE(
     +		   (XGQ(-1)+XGQ(-2)+XGQ(-2)+XGQ(-3)+XGQ(-4))*XPQ(4))+
     +               DBLE(XGQ(-4)*(XPQ(4)+XPQ(3)+XPQ(2)+XPQ(1)))
               If(IHFLA.EQ.5) then
		   QTOT =DBLE(
     +	         (XGQ(1)+XGQ(2)+XGQ(2)+XGQ(3)+XGQ(4)+XGQ(5))*XPQ(-5))+
     +               DBLE(XGQ(5)*(XPQ(-5)+XPQ(-4)+XPQ(-3)+XPQ(-2)
     + 	         +XPQ(-1)))+DBLE(
     +		   (XGQ(-1)+XGQ(-2)+XGQ(-2)+XGQ(-3)+XGQ(-4)+XGQ(-4))
     +               *XPQ(5))+DBLE(XGQ(-5)*(XPQ(5)+XPQ(4)+XPQ(3)+XPQ(2)
     +               +XPQ(1)))
               Endif
               NFLTG = -IHFLA-1
               QSUM = -QTOT*draprn()
   61          NFLTG = NFLTG + 1
               IF(NFLTG.EQ.0) GOTO 61
               IF(NFLTG.LE.0) THEN
                  NFLTP = IHFLA
               ELSE
                  NFLTP = -IHFLA
                ENDIF
               QSUM = QSUM + DBLE(XPQ(NFLTP)*XGQ(NFLTG))
               IF(QSUM.LE.0.D0.AND.NFLTG.LT.IHFLA) GOTO 61
		   IF(QSUM.LE.0) THEN
                  NFLTP = -IHFLA-1
   71             NFLTP = NFLTP + 1
                  IF(NFLTP.EQ.0) GOTO 71
                  IF(NFLTP.LE.0) THEN
                     NFLTG = IHFLA
                  ELSE
                     NFLTG = -IHFLA
                  ENDIF
                  QSUM = QSUM + DBLE(XPQ(NFLTP)*XGQ(NFLTG))
                  IF(QSUM.LE.0.D0.AND.NFLTP.LT.IHFLA) GOTO 71
		   ENDIF
c		   write(6,*) ' selected NFLTG ',NFLTG,' NFLTP ',NFLTP
      	ENDIF
            IF(NFLTG.GT.NFLAV) write(6,*) ' ELEH0 1 NFL_gam > NFLAV ',
     +      NFLTG, NFLAV
            IF(NFLTP.GT.NFLAV) write(6,*) ' ELEH0 NFL_gp > NFLAV ',
     +      NFLTP, NFLAV
            K(NIA1,2) = NFLTG
            K(NIR1,2) = -NFLTG
            K(NIA2,2) = NFLTP
            K(NF1,2) = NFLTG
            K(NF2,2) = NFLTP
         ENDIF
      ELSEIF(IRESPRO.EQ.6) THEN
C qq --> qq, q_bar q_bar --> q_bar q_bar
         SUMTF = SUMF1 + SUMF2
         IF(SUMF1/SUMTF.lt.draprn()) THEN
            ICOLORA = 2
		IF(IHFLA.LT.4) THEN
               Qtot = 0.
               Do I=1,NFLAV
                  Qtot = Qtot + XGQ(I)*XPQ(I) +  XGQ(-I)*XPQ(-I)
               End do
               QSUM = -QTOT*draprn()
               NFLT = -NFLAV-1
   80          NFLT = NFLT + 1
               IF(NFLT.EQ.0) GOTO 80
               QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(NFLT))
               IF(QSUM.LE.0.0D0) GOTO 80
		ELSE
               QTOT = (
     +	        DBLE(XGQ(IHFLA))*DBLE(XPQ(IHFLA)) + 
     +	        DBLE(XGQ(-IHFLA))*DBLE(XPQ(-IHFLA)))
               QSUM = -QTOT*draprn()
               NFLT=-IHFLA
               QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(NFLT))
               IF(QSUM.LE.0.0D0) THEN
             	NFLT=IHFLA
                  QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(NFLT))
		   ENDIF
		   IF(QSUM.Lt.0.0D0) then
		      write(6,*) ' ELEH0 problem: NFLT = ',NFLT
		   ENDIF	   
		ENDIF
c		 write(6,*) ' ELEH0: NFLT = ',NFLT
            IF(NFLT.GT.NFLAV) write(6,*) ' ELEH0 NFL > NFLAV ',NFLT,
     +      NFLAV
            K(NIA1,2) = NFLT
            K(NIR1,2) = -NFLT
            K(NIA2,2) = NFLT
            K(NF1,2) = NFLT
            K(NF2,2) = NFLT
         ELSE
            ICOLORA = 1
		IF(IHFLA.LT.4) THEN
               Qtot = 0
               Do I=1,NFLAV
                 Do J=1,NFLAV
                   Qtot = Qtot + XGQ(I)*XPQ(J) + XGQ(-I)*XPQ(-J)  
                 end do 
               end do
               NFLTG = -NFLAV-1
               QSUM = -QTOT*draprn()

   90          NFLTG = NFLTG + 1
               IF(NFLTG.EQ.0) GOTO 90
               IF(NFLTG.LE.0) THEN
                  NFLTP = -NFLAV-1
                  NFLTPC = 0
               ELSE
                  NFLTP = 0
                  NFLTPC = NFLAV
               ENDIF
  100          NFLTP = NFLTP + 1
               IF(NFLTP.GT.NFLTPC) GOTO 90
               IF(NFLTP.EQ.0) GOTO 100
               QSUM = QSUM + DBLE(XPQ(NFLTP)*XGQ(NFLTG))
               IF(QSUM.LE.0.D0.AND.NFLTP.LE.NFLTPC) GOTO 100
               IF(QSUM.LE.0.0D0) GOTO 90

		ELSE
		   QTOT =DBLE(
     +	         (XGQ(1)+XGQ(2)+XGQ(2)+XGQ(3)+XGQ(4))*XPQ(-4))+
     +               DBLE(XGQ(4)*(XPQ(-4)+XPQ(-3)+XPQ(-2)+XPQ(-1)))+
     + 	         DBLE(
     +		   (XGQ(-1)+XGQ(-2)+XGQ(-2)+XGQ(-3)+XGQ(-4))*XPQ(4))+
     +               DBLE(XGQ(-4)*(XPQ(4)+XPQ(3)+XPQ(2)+XPQ(1)))
               IF(IHFLA.EQ.5) then
		   QTOT =DBLE(
     +	         (XGQ(1)+XGQ(2)+XGQ(2)+XGQ(3)+XGQ(4)+XGQ(5))
     +               *XPQ(-5))+DBLE(XGQ(5)*(XPQ(-5)+XPQ(-4)+XPQ(-3)
     + 	         +XPQ(-2)+XPQ(-1)))+DBLE(
     +		   (XGQ(-1)+XGQ(-2)+XGQ(-2)+XGQ(-3)+XGQ(-4)+XGQ(-5))
     +               *XPQ(5))+DBLE(XGQ(-5)*(XPQ(5)+XPQ(4)+XPQ(3)+XPQ(2)
     +               +XPQ(1)))
               endif
               NFLTG = -IHFLA-1
               QSUM = -QTOT*draprn()
   91          NFLTG = NFLTG + 1
               IF(NFLTG.EQ.0) GOTO 91
               IF(NFLTG.LE.0) THEN
                  NFLTP = -IHFLA
               ELSE
                  NFLTP = IHFLA
                ENDIF
               QSUM = QSUM + DBLE(XPQ(NFLTP)*XGQ(NFLTG))
               IF(QSUM.LE.0.D0.AND.NFLTG.LT.4) GOTO 91
		   IF(QSUM.LE.0) THEN
                  NFLTP = -IHFLA-1
  101             NFLTP = NFLTP + 1
                  IF(NFLTP.EQ.0) GOTO 101
                  IF(NFLTP.LE.0) THEN
                     NFLTG = -IHFLA
                  ELSE
                     NFLTG = IHFLA
                  ENDIF
                  QSUM = QSUM + DBLE(XPQ(NFLTP)*XGQ(NFLTG))
                  IF(QSUM.LE.0.D0.AND.NFLTP.LT.IHFLA) GOTO 101
		   ENDIF
		ENDIF
c		write(6,*) ' selected NFLTG ',NFLTG,' NFLTP ',NFLTP
            IF(NFLTG.GT.NFLAV) write(6,*) ' ELEH0 2 NFL_gam > NFLAV ',
     +      NFLTG, NFLAV
            IF(NFLTP.GT.NFLAV) write(6,*) ' ELEH0 NFL_gp > NFLAV ',
     +      NFLTP, NFLAV
            K(NIA1,2) = NFLTG
            K(NIR1,2) = -NFLTG
            K(NIA2,2) = NFLTP
            K(NF1,2) = NFLTG
            K(NF2,2) = NFLTP
         ENDIF
      ELSEIF(IRESPRO.EQ.7) THEN
C qq --> qq, q_bar q_bar --> q_bar q_bar color singlet exchange
         NFLP = -NFLAV-1
         QSUM = -DBLE(XPQ2)*draprn()
  110    NFLP = NFLP + 1
         IF(NFLP.EQ.0) GOTO 110
c         write(6,*) ' q of proton ',QSUM,XPQ2,NFLP,XPQ(NFLP)
         QSUM = QSUM + DBLE(XPQ(NFLP))
         IF(QSUM.LE.0.0D0) GOTO 110
         IF(NFLP.GT.NFLAV) write(6,*) ' ELEH0 NFLP > NFLAV ',NFLP,
     +   NFLAV
         NFLG = -NFLAV-1
         QSUM = -DBLE(XGQ1)*draprn()
  120    NFLG = NFLG + 1
         IF(NFLG.EQ.0) GOTO 120
c         write(6,*) ' q of gamma ',QSUM,XGQ1,NFLG,XGQ(NFLG)
         QSUM = QSUM + DBLE(XGQ(NFLG))
         IF(QSUM.LE.0.0D0) GOTO 120
         IF(NFLG.GT.NFLAV) write(6,*) ' ELEH0 NFLG > NFLAV ',NFLG,
     +   NFLAV

         K(NIA1,2) = NFLG
         K(NIR1,2) = -NFLG
         K(NIA2,2) = NFLP
         K(NF1,2) = NFLG
         K(NF2,2) = NFLP
      ENDIF
      SCAL1 = SNGL(Q2Q)    	
      if(K(NIA1,2).EQ.21) then
         XPD1 = XGQ(0)
	   ELSE
         XPD1 = XGQ(K(NIA1,2))
	   ENDIF	
c	write(6,*) ' end of ELEH0',wt1

      IF(WT1.LE.0.D0.AND.IDEBUG.EQ.1) THEN
         write(6,*) ' ELEH0 wt1 = ',wt1,' Q2 = ',Q2,' Q2Q = ',Q2Q
         write(6,*) ' ELEH0 xgq = ',(xgq(i),i=-6,6)
         write(6,*) ' ELEH0 xpq = ',(xpq(i),i=-6,6)
      ENDIF
      IF(WT1.LT.0.D0) THEN
         write(6,*) ' ELEH0 wt1 = ',wt1,' Q2 = ',Q2,' Q2Q = ',Q2Q
         WT1 = 0.D0
      ENDIF
      RETURN
      END
