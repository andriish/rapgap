C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR A1 AS FUNCTION OF TS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSL1TS(TTS,XX,Y,XS,A1MIN,A1MAX,A1M,A1P,CFKP)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPSPC/ IPHSPC
C
      CALL HSFIVM(XX,Y,XS)
C----- A1 LIMIT FROM OMH > DELTA: A1 > A1GR
      A1GR = ( 2D0*DELTA*(PH*EPRO+PPRO*EH)
     &         - PH/XS*((XS-XX)*Y*GS + TTS - TP) )/PPRO
C
      MXS = XS*MPRO
      MXS2 = MXS*MXS
      A=((S+U-2D0*MEF2-2D0*MXS2)
     &  *(S+U-2D0*MEF2-2D0*MXS2)-4D0*T*MXS2)/16D0
      B=((U+T+S-2D0*MEF2-2D0*MXS2)
     &    *(S*TTS-U*T+(T-TTS)*(MEF2+MXS2))
     &   +T*(T-TTS)*(U-MEF2+MXS2))/8D0
      C=((TTS*(S+T)+U*T-(T+TTS)*(MEF2+MXS2))
     &  *(TTS*(S+T)+U*T-(T+TTS)*(MEF2+MXS2))
     &  -4D0*MEF2*(TTS*(U+S-2D0*MEF2-2D0*MXS2)
     &               *(U+T+S+TTS-2D0*MEF2-2D0*MXS2)
     &             +MXS2*(TTS-T)*(TTS-T)+TTS*TTS*T))/16D0
      DISK=1D0/16D0*(MEF2*(S+U-2D0*MEF2-2D0*MXS2)
     &                   *(S+U-2D0*MEF2-2D0*MXS2)
     &              +T*(MEF2+MXS2)*(S+U-MEF2-MXS2)
     &              -U*T*S + T*MXS2*(T-4D0*MEF2))
     &       *(TTS*(S+U+TTS-2D0*MEF2-2D0*MXS2)*(S+U+T-2D0*MEF2-2D0*MXS2)
     &        +MXS2*(TTS-T)*(TTS-T))
      CFKP = A
      IF (DISK.LE.0D0) THEN
        DISK = 0D0
      ENDIF
C
      IF (B.GE.0D0) THEN
        AH1 = (-B-DSQRT(DISK))/2D0/A
        AH2 = C/A/AH1
      ELSE
        AH2 = (-B+DSQRT(DISK))/2D0/A
        AH1 = C/A/AH2
      ENDIF
      A1M = DMIN1(AH1,AH2)
      A1P = DMAX1(AH1,AH2)
      A1MIN = DMAX1(A1M,A1GR)
      A1MAX = A1P
      IF (A1MIN.GE.A1MAX) THEN
        IPHSPC=1
      ELSE
        IPHSPC=0
      ENDIF
C
      END
