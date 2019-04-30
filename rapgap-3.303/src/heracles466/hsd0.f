C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSD0(GS,GT,CM1,CM2)
      COMPLEX*16 HSD0,CM1,CM2,HSSPEN
      COMPLEX*16 A,MY1,MY2,CX,CY,XX,X1X2,Y1Y2,Y1,Y2,X1,X2
      DOUBLE PRECISION GS,GT
      A = DCMPLX(-GT/GS,0D0)
      MY1 = CM1/GS
      MY2 = CM2/GS
      CX=SQRT(1D0+(MY1-MY2)*(MY1-MY2)-2D0*(MY1+MY2)+4D0*MY1*MY2/A)
      CY=SQRT(1D0+(MY1-MY2)*(MY1-MY2)-2D0*(MY1+MY2))
      XX=1D0-MY1+MY2
      X1X2=CM2/GS*(GT+CM1)/GT
      Y1Y2=CM2/GS
      IF (DREAL(XX).GT.0D0) THEN
        X1=(XX+CX)/2D0
        Y1=(XX+CY)/2D0
        X2=X1X2/X1
        Y2=Y1Y2/Y1
        ELSE
        X2=(XX-CX)/2D0
        Y2=(XX-CY)/2D0
        X1=X1X2/X2
        Y1=Y1Y2/Y2
      ENDIF
      HSD0 = 1D0/GS/GT/(X1-X2)*
     *    ( HSSPEN((1D0-X1)/(Y1-X1)) - HSSPEN(-X1/(Y1-X1))
     *     -HSSPEN((1D0-X2)/(Y2-X2)) + HSSPEN(-X2/(Y2-X2))
     *     +HSSPEN((1D0-X1)/(Y2-X1)) - HSSPEN(-X1/(Y2-X1))
     *     -HSSPEN((1D0-X2)/(Y1-X2)) + HSSPEN(-X2/(Y1-X2)) )
      END
