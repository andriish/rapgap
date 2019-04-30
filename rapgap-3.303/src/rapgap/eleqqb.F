      SUBROUTINE ELEQQB(WT1)
      IMPLICIT NONE
C
C    PHOTON GLUON ----> Q Q_BAR
C
C         P1(G)-----//////--------Q1(C)
C                  ////////
C         P2(PH)-----//////--------Q2(CBAR)
C
C
C  ( EICHLER KUNSZT ETH-PT/88-1 SOME PROBLEMS WITH AVERAGED SPIN
C  ( AND COLOUR
C  ( USE ELLIS KUNSZT CALCULATION IN NUCL PHYS. B303 (1988) 653
C  (
C

#include "rglujets.inc"
#include "rgpara1.inc"
#include "rgpara.inc"
#include "rgluco.inc"
#include "rgpart.inc"
#include "rgrapgki.inc"
      Double Precision ALPHAS,ALPHA_S,QF,SH,TH,UH,V1,V2,SUM2
      Double Precision ALPH_EM,WT1

      IF(INTER.EQ.2) THEN
         write(6,*) ' INTER = 2 not implemented for eleqqb'
         STOP
      ENDIF
      QF=DFLOAT(PYCHGE(KPA))/3.D0
      ALPHA_S=ALPHAS(DSQRT(Q2Q))
      SH = 2.D0 * AM(1)**2 + 2.D0 * DOT1(NF1,NF2)
      TH = AM(1)**2 - 2.D0 * DOT1(NIA1,NF1) - P(NIA1,5)**2
      UH = AM(2)**2 - 2.D0 * DOT1(NIA1,NF2) - P(NIA1,5)**2
c      write(6,*) ' ELEQQB ',QG2
C.....MARTIX ELEMENT
      V1 = AM(1)**2 - TH
      V2 = AM(1)**2 - UH
      SUM2 = V1/V2 + V2/V1 + 4.D0*AM(1)**2 *SH/V1/V2 -4.D0*AM(1)**4
     +*SH*SH/V1/V1/V2/V2
c      write(6,*) ' old sum2 : ',sum2
c      write(6,*) ' sss ',sss
cccccccccccccc         SUM2 = SUM2/QG2**2

      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(Q2)

      SUM2 = SUM2*16.D0*PI*PI*ALPHA_S*ALPH_EM*QF*QF
      WT1=SUM2
      RETURN
      END
