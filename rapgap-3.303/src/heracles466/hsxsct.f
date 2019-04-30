C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   XSCUT FROM OMEGA-MIN(XS)=DELTA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSXSCT(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      XS1 = X
      XS2 = 1D0
      DO 100, N = 1,70,1
        XS3 = (XS1+XS2)/2D0
        OM3 = HSOMIN(X,Y,XS3)
        IF (OM3.LT.DELTA) THEN
          XS1 = XS3
        ELSE
          XS2 = XS3
        ENDIF
100   CONTINUE
      HSXSCT = XS3
      END
