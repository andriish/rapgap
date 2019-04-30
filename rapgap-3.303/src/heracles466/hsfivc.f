C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   CAPITAL INVARIANTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIVC(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
C
      GS = SP-MEI2-MPRO2
      TP = -X*Y*GS
      GU = -(1D0-Y)*GS
      UP = GU+MPRO2+MEF2
      END
