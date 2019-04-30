C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SOFT BREMSSTRAHLUNG FOR ELECTRON-QUARK SCATTERING
C
      FUNCTION HSBRNC(X,Y,LL,LQF)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSBRNC,HSSPEN
     1       ,BEIEF,BQIQF,BEIQI,BEFQF,BEFQI,BEIQF
     2       ,ZEIEF,ZQIQF,ZEIQI,ZEFQF,ZEFQI,ZEIQF
      DIMENSION CQFL(2)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSISGM/ TCUTQ,TCUTQS
      DATA CQFL /0.6666666666666666D0, -0.3333333333333333D0/
C
      CQF=CQFL(LQF)
      GSP=SP-MEI2-MPRO2
      Q2=GSP*X*Y
      S=X*SP
      U=-(S-Q2)
      DELTA2=DELTA*DELTA
      EEI=EELE
      EEI2=EEI*EEI
      EQI=SP/4D0/EEI*X
      EQI2=EQI*EQI
      EEF=EEI*(1D0-Y)+Q2/4D0/EEI
      EEF2=EEF*EEF
      EQF=EEI+EQI-EEF
      EQF2=EQF*EQF
C
      HSBRNC=0D0
      IF (LPAR(12).EQ.1) THEN
        ZEIEF=DCMPLX(1D0-4D0*EEI*EEF/Q2,0D0)
        BEIEF=-HSSPEN(ZEIEF)
        HSBRNC=HSBRNC
     *          -2D0*DLOG(MEI2/Q2)*(1D0+DLOG(4D0*DELTA2/Q2) )
     *          -DLOG(DELTA2/EEI2)-DLOG(DELTA2/EEF2)
     *          +2D0*DREAL(BEIEF)-2D0*PI*PI/3D0
     *          -DLOG(EEI2/EEF2)*DLOG(EEI2/EEF2)/4D0
     *          +DLOG(Q2/4D0/EEI/EEF)
     *             *(DLOG(4D0*EEI2/MEI2)+DLOG(4D0*EEF2/MEI2)
     *               +DLOG(Q2/4D0/EEI/EEF)                )
      ENDIF
      IF (LPAR(13).EQ.1) THEN
        TCQ2=TCUTQ*TCUTQ
        TCQS2=TCUTQS*TCUTQS
        ZQIQF=DCMPLX(1D0-4D0*EQI*EQF/Q2,0D0)
        BQIQF=-HSSPEN(ZQIQF)
        HSBRNC=HSBRNC+CQF*CQF*(
     *           -DLOG(Q2/MQI2)
     *           +2D0*DREAL(BQIQF)
     *           +9D0/2D0-4D0*PI*PI/3D0
     *           -DLOG(EQI2/EQF2)*DLOG(EQI2/EQF2)/4D0
     *           -DLOG(4D0*EQI*EQF/Q2)*DLOG(4D0*EQI*EQF2/Q2)
     *           -(3D0/2D0+DLOG(DELTA2/EQF2))
     *                 *DLOG((TCQS2*EQF2+MQF2)/Q2)
     *           -(3D0/2D0+DLOG(DELTA2/EQI2))
     *                 *DLOG((TCQ2*EQI2+MQI2)/Q2))
      ENDIF
      IF (LPAR(14).EQ.1) THEN
        ZEIQI=DCMPLX(1D0-4D0*EEI*EQI/S,0D0)
        BEIQI=-HSSPEN(ZEIQI)
        ZEFQF=DCMPLX(1D0-4D0*EEF*EQF/S,0D0)
        BEFQF=-HSSPEN(ZEFQF)
        ZEFQI=DCMPLX(1D0-4D0*EEF*EQI/(-U),0D0)
        BEFQI=-HSSPEN(ZEFQI)
        ZEIQF=DCMPLX(1D0-4D0*EEI*EQF/(-U),0D0)
        BEIQF=-HSSPEN(ZEIQF)
        HSBRNC=HSBRNC+2D0*LL*CQF*(
     *           +2D0*DLOG(4D0*DELTA2/Q2)*DLOG((-U)/S)
     *           +DREAL(-BEIQI-BEFQF+BEFQI+BEIQF)      )
      ENDIF
      RETURN
      END
