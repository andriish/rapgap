      SUBROUTINE ELEQQL(WT1)
C
C    PHOTON GLUON ----> Q Q_BAR
C
C         P1(G)-----//////--------Q1(q)
C                  ////////
C         P2(PH)-----//////--------Q2(q_bar)
C
C
C
C
C for light quarks without masses
c changed for including Q2 of photon
C
      IMPLICIT NONE
#include "rglujets.inc"
#include "rgpara1.inc"
#include "rgpara.inc"
#include "rgluco.inc"
#include "rgrapgki.inc"
      Double Precision WMAX
	Integer IMIX
      COMMON /OALPHAS/ WMAX,IMIX
      Double Precision P2E(4),Q1E(4),Q2E(4)
	Double Precision ALPHA_S,ALPHAS,ALPH_EM,SH,TH,UH,QF2
	Double Precision SUMT,SUML,EPSILON,YX,WT1,DOT
	Integer I 
      IF(INTER.EQ.2) THEN
         write(6,*) ' interaction INTER = 2 not implemented '
         STOP
      ENDIF
      DO 10  I=1,4
         P2E(I) = P(NIA2,I)
         Q1E(I) = P(NF1,I)
         Q2E(I) = P(NF2,I)
   10 CONTINUE
   20 CONTINUE
C.....MARTIX ELEMENT
      ALPHA_S=ALPHAS(DSQRT(Q2Q))
      SH =  2.D0 * DOT(Q1E,Q2E)
      TH = - 2.D0 * DOT(P2E,Q2E)
      UH = - 2.D0 * DOT(P2E,Q1E)
c sum of charges for u d s (2/3)**2 + (1/3)**2 + (1/3)**2 = 6/9 = 2/3
      QF2=2.D0/3.D0
      IF(IMIX.EQ.1) QF2 = DFLOAT(PYCHGE(KPA))**2/9.D0
c for transverse photons
      SUMT = 16.D0*PI*PI*(UH/TH + TH/UH
     +        - 2.D0*SH*Q2/UH/TH
     +        + 4.D0*SH*Q2/(SH+Q2)**2)
c for longitudinal photons
      SUML = -16.D0*8.D0*PI*PI*Q2*(Q2+TH+UH)/(SH+Q2)**2
      IF(SUMT.LE.0.0D0.OR.SUML.LT.0.0D0) THEN
         write(6,*) ' eleqql ',TH,SH,UH
         write(6,*) ' eleqql ',SUMT,SUML
      endif

c epsilon
      YX = DBLE(YY)
      EPSILON = (1.D0 - YX)/(1.D0 - YX + YX**2 /2)
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)

      WT1 = ALPHA_S*ALPH_EM*QF2*(SUMT + EPSILON * SUML)
      IF(WT1.LE.0.D0) THEN
         write(6,*) ' eleqql ',WT1,EPSILON
         write(6,*) ' eleqql ',TH,SH,UH
         write(6,*) ' eleqql ',SUMT,SUML
      endif

      RETURN
      END
