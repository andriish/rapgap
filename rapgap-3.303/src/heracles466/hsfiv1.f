C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INVARIANTS CALCULATED IN KP - CMS - SYSTEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIV1(X,Y,XS,A1,TTS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      CALL HSFIVM(X,Y,XS)
      DKP = A1/2D0
      DKPS = (TTS-T+A1)/2D0
      DKQS =(XS-X)*Y*GS/2D0
      DKQ = DKQS +(TTS-T)/2D0
      SS = S - 2D0*(DKP+DKQ)
      US = U + 2D0*(DKPS-DKQ)
      GX = -TTS/XS
      END
