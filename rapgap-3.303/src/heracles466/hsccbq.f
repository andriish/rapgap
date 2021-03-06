C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCBQ(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSSPEN,ZEQI,ZEQF,ZQIQF,CEQI,CEQF,CQIQF
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      HSCCBQ=0D0
      CQU=2D0/3D0
      GSP=SP-MEI2-MPRO2
      S=X*SP
      T=-GSP*X*Y
      U=-X*SP-T
      EEI=EELE
      EQI=X*SP/4D0/EEI
      ENU=EEI*(1D0-Y)-T/4D0/EEI
      EQF=EEI+EQI-ENU
C
      DLMEEE=DLOG(ME2/4D0/EEI/EEI)
      DLMUEI=DLOG(MU2/4D0/EQI/EQI)
      DLMDEF=DLOG(MD2/4D0/EQF/EQF)
      DLMWME=DLOG(MW2/ME2)
      DLMWMU=DLOG(MW2/MU2)
      DLMWMD=DLOG(MW2/MD2)
      DLQ2MW=DLOG(-T/MW2)
C
C---LEPTONIC PART
      IF (LPAR(12).EQ.1) THEN
        CEQF=DCMPLX(1D0+4D0*EEI*EQF/U,0D0)
        ZEQF=HSSPEN(CEQF)
        REQF=-DREAL(ZEQF)-DLMEEE*DLMEEE/4D0-DLMDEF*DLMDEF/4D0-PI*PI/3D0
        HSCCBQ=HSCCBQ+ALP2PI
     .  *(2D0*REQF-DLOG(DELTA*DELTA/EQF/EQF)-DLOG(DELTA*DELTA/EEI/EEI)
     .    +0.5D0*DLMWME*(DLMWME+3D0)+0.5D0*DLMWMD*(DLMWMD+3D0)
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-ME2/U)+DLOG(-MD2/U))
     .    -DLQ2MW*(DLQ2MW-2D0*DLOG(T/U)))
      ENDIF
C---QUARKONIC PART
      IF (LPAR(13).EQ.1) THEN
        CQIQF=DCMPLX(1D0+4D0*EQI*EQF/T,0D0)
        ZQIQF=HSSPEN(CQIQF)
        RQIQF=-DREAL(ZQIQF)-DLMUEI*DLMUEI/4D0-DLMDEF*DLMDEF/4D0
     .        -PI*PI/3D0
        HSCCBQ=HSCCBQ+ALP2PI*CQU*CQU
     .  *(2D0*RQIQF+0.5D0*DLMWMU*(DLMWMU+3D0)+0.5D0*DLMWMD*(DLMWMD+3D0)
     .    -DLOG(DELTA*DELTA/EQI/EQI)-DLOG(DELTA*DELTA/EQF/EQF)
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MU2/T)+DLOG(-MD2/T))
     .    -DLQ2MW*(DLQ2MW-3D0))
      ENDIF
C---LEPTON QUARK INTERFERENCE
      IF (LPAR(14).EQ.1) THEN
        CEQF=DCMPLX(1D0+4D0*EEI*EQF/U,0D0)
        ZEQF=HSSPEN(CEQF)
        REQF=-DREAL(ZEQF)-DLMEEE*DLMEEE/4D0-DLMDEF*DLMDEF/4D0-PI*PI/3D0
        CQIQF=DCMPLX(1D0+4D0*EQI*EQF/T,0D0)
        ZQIQF=HSSPEN(CQIQF)
        RQIQF=-DREAL(ZQIQF)-DLMUEI*DLMUEI/4D0-DLMDEF*DLMDEF/4D0
     .        -PI*PI/3D0
        CEQI=DCMPLX(1D0-4D0*EEI*EQI/S,0D0)
        ZEQI=HSSPEN(CEQI)
        REQI=-DREAL(ZEQI)-DLMEEE*DLMEEE/4D0-DLMUEI*DLMUEI/4D0-PI*PI/3D0
        HSCCBQ=HSCCBQ+ALP2PI*(-CQU)
     .  *(-2D0*(REQI-RQIQF-REQF)-2D0*DLOG(DELTA*DELTA/EQF/EQF)
     .    +DLMWMD*(DLMWMD+3D0)
     .    +DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(ME2/S)+DLOG(MU2/S))
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MU2/T)+DLOG(-MD2/T))
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-ME2/U)+DLOG(-MD2/U))
     .    -DLQ2MW*(DLQ2MW-3D0+2D0*DLOG(-U/S)))
      ENDIF
      END
