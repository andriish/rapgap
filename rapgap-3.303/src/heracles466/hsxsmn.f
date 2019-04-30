C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   XSMIN FROM OMEGA-MAX(XS)=DELTA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSXSMN(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      XS1 = X
      XS2 = 1D0
      DO 90, N = 1,70,1
        XS3 = (XS1+XS2)/2D0
        OM3 = HSOMAX(X,Y,XS3)
        IF (OM3.LT.DELTA) THEN
          XS1 = XS3
        ELSE
          XS2 = XS3
        ENDIF
90    CONTINUE
      HSXSMN = XS3

C...CHECK ON ES > ME
      IF (X.LT.1D-6) THEN
      BXS=GU*(TP+MEI2)/(GU*GU-4D0*MEI2*MPRO2)
      DXS=1D0-(GU*GU-4D0*MEI2*MPRO2)/GU/GU
     *        *(TP*(TP-2D0*MEI2)-7D0*MEI2*MEI2)/(TP+MEI2)/(TP+MEI2)
      IF (DXS.LT.0D0) RETURN
      XESMIN=-BXS*(1D0+DSQRT(DXS))
      XESMAX=-BXS*(1D0-DSQRT(DXS))
      HSXSMN=DMAX1(HSXSMN,XESMIN,XESMAX)
      ENDIF
      END
