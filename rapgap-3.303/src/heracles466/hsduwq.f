C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSDUWQ(T)
C---IR-FINITE PART
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSDUWQ,HSCLN,HSCLM1,HSCLM2,HSCLM4,CMW2,CMZ2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      CQU=2D0/3D0
      HSDUWQ=DCMPLX(0D0,0D0)
C
C---PHOTONIC PART
      IF (LPAR(12).EQ.1) THEN
        HSDUWQ=HSDUWQ+ALP4PI*(-HSCLN(CMZ2/MD2)-9D0/2D0
     *                        +3D0*HSCLM4(T,CMW2,(0D0,0D0),MD2))
      ENDIF
      IF (LPAR(14).EQ.1) THEN
        AMBU=HSCLM1(T,MU2)+DLOG(-T/MU2)*DLOG(-T/MU2)
        AMBD=HSCLM1(T,MD2)+DLOG(-T/MD2)*DLOG(-T/MD2)
        HSDUWQ=HSDUWQ+ALP4PI*CQU*(
     *      -(3D0*DLOG(MD/MU)+.5D0*AMBU+.5D0*AMBD)
     *      -SW2/CW2*HSCLM2(T,CMZ2)
     *      +3D0*HSCLM4(T,CMW2,(0D0,0D0),MU2)
     *      -3D0*HSCLM4(T,CMW2,(0D0,0D0),MD2)
     *      +HSCLN(CMZ2/MD2)+9D0/2D0)
      ENDIF
      IF (LPAR(13).EQ.1) THEN
        AMBU=HSCLM1(T,MU2)+DLOG(-T/MU2)*DLOG(-T/MU2)
        AMBD=HSCLM1(T,MD2)+DLOG(-T/MD2)*DLOG(-T/MD2)
        HSDUWQ=HSDUWQ+ALP4PI*CQU*CQU*(
     *             3D0*DLOG(MD/MU)+.5D0*AMBU+.5D0*AMBD
     *             +SW2/CW2*HSCLM2(T,CMZ2))
      ENDIF
C---WEAK PART
      IF (LPAR(15).EQ.1) THEN
        HSDUWQ=HSDUWQ+ALP4PI*(
     *      +(2D0*SW2-1D0)/4D0/SW2/CW2*HSCLM2(T,CMZ2)
     *      +3D0*CW2/SW2*HSCLM4(T,CMW2,CMZ2,0D0)
     *      +3D0/SW2+(1D0/2D0/SW2-3D0*CW2/SW2/SW2)*DLOG(1D0/CW2))
      ENDIF
      END
