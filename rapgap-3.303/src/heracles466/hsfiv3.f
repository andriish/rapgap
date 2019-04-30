C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INVARIANTS CALCULATED IN KQ - CMS - SYSTEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIV3(X,Y,XS,A1,A3)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGIKP/ GS,GU,GX,TP,UP

      CALL HSFIVM(X,Y,XS)
      DKP = A1/2D0
      DKQ = A3/2D0
      DKQS =(XS-X)*Y*GS/2D0
      DKPS = DKP+DKQ-DKQS
      TS = T + A3 - 2D0*DKQS
      SS = S - A1 - A3
      US = U + A1 - 2D0*DKQS
      GX = -TS/XS
      END
