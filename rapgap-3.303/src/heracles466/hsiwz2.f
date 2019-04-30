C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       FUNCTION HSIWZ2(GS,GT,CM1,CM2)
       COMPLEX*16 HSIWZ2,HSD0,HSCMW,HSCMWZ,CM1,CM2,HSCLN,HSFONE
       COMPLEX*16 A,MY1,MY2,C,Y1,Y2,COEFF1,COEFF2,COEFF3
       DOUBLE PRECISION GS,GT,RM1,RM2
       A = DCMPLX(-GT/GS,0D0)
       MY1 = CM1/GS
       MY2 = CM2/GS
       RM1=DSQRT(DREAL(CM1))
       RM2=DSQRT(DREAL(CM2))
       C=SQRT(1D0+(MY1-MY2)*(MY1-MY2)-2D0*(MY1+MY2))
       Y1=(1D0-MY1+MY2+C)/2D0
       Y2=(1D0-MY1+MY2-C)/2D0
       COEFF1=0.25D0/(GS+GT)/(GS+GT)*
     *      ( GS*(GS*GT+2D0*CM1*CM2)
     *       -GT*(GS-CM1)*(GS-CM1) - GT*(GS-CM2)*(GS-CM2)
     *       -2D0*GT*(GS+GT)*(GT+CM1+CM2)                )
       COEFF2=-GT*(GS+2D0*GT+CM1+CM2)/4D0/(GS+GT)/(GS+GT)
       COEFF3=(GS*(GS-CM1-CM2)+2D0*GT*(GS+GT))/4D0/(GS+GT)/(GS+GT)
       HSIWZ2=COEFF1*HSD0(GS,GT,CM1,CM2)
     *     +COEFF2*(HSCMW(GT,CM1)+HSCMW(GT,CM2))
     *     +COEFF3*2D0*HSCMWZ(GS,CM1,CM2)
     *     -0.5D0/(GS+GT)*( 1D0 - HSCLN(-GT/CM1)
     *                      -CM2/(CM1-CM2)*HSCLN(CM2/CM1)
     *                      -HSFONE(GS,RM1,RM2)             )
       END
