C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   XSMIN FROM OMEGA-MAX(XS)=DELTA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSXMAX(Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
C
      XS=1D0
      X1=XMIN
      X2=1D0
      DO 90 N=1,70,1
        X3=(X1+X2)/2D0
        CALL HSFIVC(X3,Y)
        OM3=HSOMAX(X3,Y,XS)
        IF (OM3.LT.DELTA) THEN
          X2=X3
        ELSE
          X1=X3
        ENDIF
90    CONTINUE
      HSXMAX=X3

C...CHECK ON ES > ME
C     IF (X3.LT.1D-6) THEN
C     BXS=GU*(TP+MEI2)/(GU*GU-4D0*MEI2*MPRO2)
C     DXS=1D0-(GU*GU-4D0*MEI2*MPRO2)/GU/GU
C    *        *(TP*(TP-2D0*MEI2)-7D0*MEI2*MEI2)/(TP+MEI2)/(TP+MEI2)
C     IF (DXS.LT.0D0) RETURN
C     XEMIN=-BXS*(1D0+DSQRT(DXS))
C     XEMAX=-BXS*(1D0-DSQRT(DXS))
C     HSXMAX=DMIN1(HSXMAX,XEMIN,XEMAX)
C     ENDIF

      END
