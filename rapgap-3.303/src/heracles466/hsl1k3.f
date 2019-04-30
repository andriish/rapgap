C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR A1 AS FUNCTION OF A3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSL1K3(A3,XX,Y,XS,A1MIN,A1MAX,A1M,A1P,CFKP)
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
      IPHSPC=0
      CALL HSFIVM(XX,Y,XS)

C---A1 LIMIT FROM OMH > DELTA: A1 > A1GR
      A1GR = (2D0*DELTA*(PH*EQH+PQH*EH) - PH*A3)/PQH

      A = ((S+U-2D0*MEF2-2D0*MQF2)*(S+U-2D0*MEF2-2D0*MQF2) - 4D0*T*MQF2)
     &    / 16D0
      B =( (U+S+T-A3-2D0*MEF2-2D0*MQF2)
     &   *(2D0*T*MQF2 - (S+U-2D0*MEF2-2D0*MQF2)
     &                 *(S-MEF2-MQF2) )
     &    + T*A3*(S-U) )/8D0
      C=( ((U+S+T-A3-2D0*MEF2-2D0*MQF2)*(S-MEF2-MQF2)-A3*T)
     &   *((U+S+T-A3-2D0*MEF2-2D0*MQF2)*(S-MEF2-MQF2)-A3*T)
     &  -4D0*MEF2*( MQF2*(S+U+T-2D0*MEF2-2D0*MQF2)
     &                  *(S+U+T-2D0*MEF2-2D0*MQF2)
     &       -A3*(S+U-A3-2D0*MEF2-MQF2)*(S+U+T-2D0*MEF2-MQF2)
     &       -A3*MQF2*(T-MQF2)) )/16D0
      DISK=( MEF2*(U+S-2D0*MEF2-2D0*MQF2)
     &           *(U+S-2D0*MEF2-2D0*MQF2)
     &         +T*(U+S-MEF2-MQF2)*(MEF2+MQF2)
     &         -U*T*S + T*MQF2*(T-4D0*MEF2) )
     &   * ( MQF2*(S+U+T-2D0*MEF2-2D0*MQF2)
     &           *(S+U+T-2D0*MEF2-2D0*MQF2)
     &        -A3*(S+U+T-2D0*MEF2-2D0*MQF2)*(S+U-A3-2D0*MEF2)
     &        + A3*A3*MQF2 )/16D0

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
      ENDIF
      END
