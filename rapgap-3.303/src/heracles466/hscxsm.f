C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   XSMIN FROM OMEGA-MAX(XS)=DELTA (FOR CC)
C   FROM CXSMIN (epcctot.f, 28/08/98)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCXSM(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      XS1=X+(MEI2+MQI2+MQF2)/SP/Y
      XSMM=XS1
      XS2=1D0
      DO 90, N=1,70,1
        XS3=(XS1+XS2)/2D0
        OM3=HSCKMX(X,Y,XS3)
        IF (OM3.LT.DELTA) THEN
          XS1=XS3
        ELSE
          XS2=XS3
        ENDIF
90    CONTINUE
      HSCXSM=DMAX1(XS3,XSMM)
      END
