CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      FUNCTION HSSRWW(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSRWW(Q2): BEITRAEGE ZUR Z-BOSONSELBSTENERGIE (RENORMIERT)
C              Q2 = Q**2 VIEREIMPULSQUADRAT
C-----------------------------------------------------------------------
C       21.10.83
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSSRWW,CJFINW,CAFINZ,CAFINW
     *          ,HSSFZZ,HSSFWW
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
        SQ=ALP4PI*(DLOG(MU2/MD2)+DLOG(MC2/MS2)+DLOG(MT2/MB2))
        SP=ALP4PI*DLOG(MT2/MB2)*(MT2-MB2)/4./SW2/MW2
        CAFINZ=HSSFZZ(MZ2)
        CAFINW=HSSFWW(MW2)
        RDMW2=DREAL(CAFINW)
        RDMZ2=DREAL(CAFINZ)
        MRENKW=CW2/SW2*(RDMZ2/MZ2 - RDMW2/MW2) + ALP2PI/3D0
     *         + SQ/3./SW2 + SP*CW2/SW2
        CJFINW=DCMPLX(RDMW2,0D0)
        HSSRWW=(HSSFWW(Q2)-CJFINW+(Q2-MW2)*MRENKW)
        IF (LPAR(7).EQ.2) THEN
          DDALPP=HSDSGQ(MZ2)*(Q2-MW2)
          HSSRWW=HSSRWW+DCMPLX(DDALPP,0D0)
        ENDIF
      RETURN
      END