CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE HSSETF(T)
C---FORM FACTORS
C   RECALCULATE EFFECTIVE COUPLING CONSTANTS INCLUDING SELF ENERGIES
C   RUNNING ALPHA, EFFECTIVE SW2, FORM FACTOR KAPPA (FROM PI_Z)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSSRGG,HSSRGZ,HSSRZZ,CG,CM,CZ,HSFHFB
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSDELR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMC1/ VAFI1(2,3,2)
     &               ,AFIJ1(3,2,2),BFIJ1(3,2,2),FLIND1(2,3,2,2)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSKNCC/ SXNRCC,SX1NCC
      COMMON /HSFRFF/ ALPFFQ,AKAPPA,GMUFFQ,SWEFF2
C
C  SELF ENERGIES
C
      IF (LPAR(7).GE.1) THEN
        CG=HSSRGG(T)
        PIGGG=DREAL(CG)/T
        ELSE
        PIGGG=0D0
      ENDIF
      ALPFFQ=1D0/(1D0+PIGGG)
      EFFQ=SQRT(ALPFFQ)
      IF (LPAR(8).EQ.1) THEN
        CM=HSSRGZ(T)
        SM=DREAL(CM)/T
        ELSE
        SM=0D0
      ENDIF
      AKAPPA=1D0-CW/SW*SM/(1D0+PIGGG)
      SWEFF2=SW2*AKAPPA
      IF (LPAR(9).EQ.1) THEN
        CZ=HSSRZZ(T)
        PIGZZ=DREAL(CZ)/(T-MZ2)
        ELSE
        PIGZZ=0D0
      ENDIF
      GMUFFQ=1D0/(1D0+PIGZZ)
      BFFQ=SQRT(GMUFFQ)
      IF (LPAR(3).GE.1) BFFQ=BFFQ*BTOP4
C
C---DEFINE FERMION GAUGE BOSON COUPLING CONSTANTS
      IF (LPAR(4).EQ.1) THEN
        B0=1D0/4D0/CW/SW
        B=B0
        ELSE
C---NORMALIZED TO G-MU
        B0=MZ/SQRT(AGF0)/4D0
        B=B0
        IF (LPAR(9).GE.1) B=B0*SQRT(1D0-DELTAR)
      ENDIF
      IF (LPAR(2).EQ.1.AND.LPAR(9).GE.1) B=B*BFFQ
      RHO1=B/B0

      VAFI1(2,1,1)=0D0
      VAFI1(2,2,1)=0D0
      VAFI1(2,3,1)=0D0
      VAFI1(2,1,2)=-B
      VAFI1(2,2,2)=B
      VAFI1(2,3,2)=-B
      VAFI1(1,1,1)=EFFQ
      VAFI1(1,2,1)=-2D0/3D0*EFFQ
      VAFI1(1,3,1)=1D0/3D0*EFFQ
      VAFI1(1,1,2)=B*(4D0*SWEFF2-1D0)
      VAFI1(1,2,2)=B*(1D0-8D0*SWEFF2/3D0)
      VAFI1(1,3,2)=B*(4D0*SWEFF2/3D0-1D0)
C
C..VERTEX CORRECTIONS (NO QED PARTS)
      LPRK11=LPAR(11)
      LPAR(11)=0
C..ELECTRON VERTEX
      IF (LPAR(12).EQ.1) THEN
        LF=1
        DO 1 IVA=1,2
        VAFI1(IVA,LF,1)=VAFI1(IVA,LF,1)
     *                   +EFFQ*DREAL(HSFHFB(T,IVA,LF,1,ME2))
        VAFI1(IVA,LF,2)=VAFI1(IVA,LF,2)
     *                   +RHO1*DREAL(HSFHFB(T,IVA,LF,2,ME2))
    1   CONTINUE
      ENDIF
      LPAR(11)=LPRK11
C..QUARK VERTEX
      IF (LPAR(13).EQ.1) THEN
        DO 2 LF=2,3
        DO 2 IVA=1,2
        VAFI1(IVA,LF,1)=VAFI1(IVA,LF,1)
     *                   +EFFQ*DREAL(HSFHFB(T,IVA,LF,1,MQI2))
        VAFI1(IVA,LF,2)=VAFI1(IVA,LF,2)
     *                   +RHO1*DREAL(HSFHFB(T,IVA,LF,2,MQI2))
    2   CONTINUE
      ENDIF
C
      IGAMMA = 1
      IZ = 2
      IEL = 1
      IFU = 2
      IFD = 3
      INDV = 1
      INDA = 2
C
      DO 3 IF=IEL,IFD
        DO 3 IB1=IGAMMA,IZ
          DO 3 IB2=IGAMMA,IZ
          FLIND1(INDV,IF,IB1,IB2)=
     *     2D0*(VAFI1(INDV,IF,IB1)*VAFI1(INDV,IF,IB2)
     *         +VAFI1(INDA,IF,IB1)*VAFI1(INDA,IF,IB2))
    3 CONTINUE
      DO 4 IF=IEL,IFD
        DO 4 IB1=IGAMMA,IZ
          DO 4 IB2=IGAMMA,IZ
          FLIND1(INDA,IF,IB1,IB2)=
     *     2D0*(VAFI1(INDV,IF,IB1)*VAFI1(INDA,IF,IB2)
     *         +VAFI1(INDA,IF,IB1)*VAFI1(INDV,IF,IB2))
    4 CONTINUE
C
      DO 5 IVB1 = IGAMMA, IZ
       DO 5 IVB2 = IGAMMA, IZ
        DO 5 IFERM = IFU, IFD
        AFIJ1(IFERM,IVB1,IVB2)=FLIND1(INDV,IFERM,IVB1,IVB2)
     *  *(FLIND1(INDV,IEL,IVB1,IVB2)-POLARI*FLIND1(INDA,IEL,IVB1,IVB2))
        BFIJ1(IFERM,IVB1,IVB2)=FLIND1(INDA,IFERM,IVB1,IVB2)
     *  *(FLIND1(INDA,IEL,IVB1,IVB2)-POLARI*FLIND1(INDV,IEL,IVB1,IVB2))
    5 CONTINUE
C
      RETURN
      END