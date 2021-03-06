C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR TS FOR KPS - PART
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSLTS2(A2,XX,Y,XS,TSMIN,TSMAX,TSM,TSP,CFKPS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPSPC/ IPHSPC
C
      CALL HSFIVM(XX,Y,XS)
C---TS LIMIT FROM OMH > DELTA
      TGR=2D0*XS*DELTA*(PH*EPRO+PPRO*EH)
     &   -XS*PPRO*(TP+A2)-PH*((XS-XX)*Y*GS-TP)
C---CUT IN HADRONIC MASS
      IF (ICUT.GT.1.AND.WMIN.GT.MPRO) THEN
        W2CUT=WMIN*WMIN
        TWCUT=-XS*(W2CUT-MPRO2)/(1D0-XS)
      ELSE
        TWCUT=0D0
      ENDIF
C
      MXS=XS*MPRO
      MXS2=MXS*MXS
      A=((U+T-3D0*MEF2-MXS2)*(U+T-3D0*MEF2-MXS2)-4D0*MEF2*S)/16D0
      B=(-2D0*MEF2*(U+S-2D0*MEF2-2D0*MXS2)*(U+S-2D0*MEF2-2D0*MXS2)
     &   +(U+S+T-2D0*MEF2-2D0*MXS2)
     &     *((A2-T)*(MEF2+MXS2)-2D0*T*MEF2-U*A2)
     &   +T*(S*(U+T+A2)-(MEF2-MXS2)*(MEF2-MXS2)
     &       -A2*(MEF2-MXS2)+2D0*T*MEF2))/8D0
      C=((A2*(U+S-2D0*MEF2-2D0*MXS2)+T*(S-MEF2-MXS2))
     &  *(A2*(U+S-2D0*MEF2-2D0*MXS2)+T*(S-MEF2-MXS2))
     &  -4D0*T*T*MXS2*(A2+MEF2) -4D0*T*A2*A2*MXS2)/16D0
      DISK=1D0/16D0*(MEF2*(S+U+T-2D0*MEF2-2D0*MXS2)
     &                   *(S+U+T-2D0*MEF2-2D0*MXS2)
     &               +A2*(S+U+T-2D0*MEF2-2D0*MXS2)*(A2+U+T-MEF2-MXS2)
     &               +A2*A2*MXS2)
     &       *(MEF2*(S+U-2D0*MEF2-2D0*MXS2)
     &             *(S+U-2D0*MEF2-2D0*MXS2)
     &         +T*(MEF2+MXS2)*(S+U-MEF2-MXS2)
     &         -U*T*S+T*MXS2*(T-4D0*MEF2))
      CFKPS=A
      IF (DISK.LE.0D0) THEN
        DISK=0D0
      ENDIF
C
      IF (B.GE.0D0) THEN
        T1=(-B-DSQRT(DISK))/2D0/A
        T2=C/A/T1
      ELSE
        T2=(-B+DSQRT(DISK))/2D0/A
        T1=C/A/T2
      ENDIF
      TSM=DMIN1(T1,T2)
      TSP=DMAX1(T1,T2)
C
      TSMIN=TSM
      TSMAX=TSP
C
      IF (PH.GT.(XS*PPRO)) THEN
        TGR=TGR/(PH-XS*PPRO)
        TSMIN=DMAX1(TSMIN,TGR)
      ENDIF
      IF (PH.LT.(XS*PPRO)) THEN
        TGR=TGR/(PH-XS*PPRO)
        TSMAX=DMIN1(TGR,TSMAX)
      ENDIF
      TSMAX=DMIN1(TSMAX,TWCUT)
C
      IF (TSMIN.GE.TSMAX) THEN
        IPHSPC=1
      ELSE
        IPHSPC=0
      ENDIF
      END
