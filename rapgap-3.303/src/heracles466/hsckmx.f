C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   MAXIMAL PHOTON ENERGY AS FUNCTION OF XS (FOR CC)
C   FROM CK0MAX (epcctot.f, 28/08/98)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCKMX(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
C
      S=XS*SP
      T=-X*Y*SP
      U=-XS*(1D0-Y)*SP
      SIGMA=(S-MEI2-MQI2)/2D0
      TAU=-(T-MEI2)/2D0
      QRHO=-(U-MQI2)/2D0
      QH1=-(SIGMA*PELE+EELE*DSQRT(SIGMA*SIGMA-MEI2*MQI2))
      QPQH=-(SIGMA*SIGMA-EELE*EELE*MQI2)/QH1
      QEQH=DSQRT(QPQH*QPQH+MQI2)
      EESH=(QPQH*TAU+PELE*QRHO)/(PELE*QEQH+QPQH*EELE)
      PPSH=EESH
      ECOS=(QRHO*EELE-TAU*QEQH)/PPSH/(PELE*QEQH+QPQH*EELE)
      ESIN=DSQRT(1D0-ECOS*ECOS)
      SMF=MEI2+MQI2+MQF2
      STU=S+T+U-SMF
      STUM=STU+MQF2
      CK0MAX=STU/2D0/STUM*(EELE+QEQH-EESH
     *  +DSQRT(PPSH*PPSH*ESIN*ESIN
     *         +(PELE-QPQH-PPSH*ECOS)*(PELE-QPQH-PPSH*ECOS)))
      HSCKMX=CK0MAX
      END