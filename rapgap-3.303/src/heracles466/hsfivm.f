C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   MANDELSTAM INVARIANTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIVM(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      S = XS*GS + XS*XS*MPRO2 + MEI2
      T = -X*Y*GS
      U = -XS*(1D0-Y)*GS + XS*XS*MPRO2 + MEF2
      END
