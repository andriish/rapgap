C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE HSSAB1(X,Y,LL,POL,LQ,RA1,RA3,RB1,RB3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     1-LOOP VIRTUAL CORRECTIONS TO DEEP INELASTIC ELECTRON QUARK
C     SCATTERING
C     LL = +1,-1: POSITRONS, ELECTRONS
C     POL  = +1,-1: RIGHT-, LEFTHANDED LEPTONS
C
C--->
C--->  NEW VERSION 07.02.1991
C--->  FORM FACTOR APPROACH
C--->
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 CMZ2,CMW2
     *       ,HSBRNC,BSQ,HSBCGA,HSBXCV,HSBXCA,HSBXI0,HSBXI5
     *       ,CVEG,HSFHFB
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMC1/ VAFI1(2,3,2)
     &               ,AFIJ1(3,2,2),BFIJ1(3,2,2),FLIND1(2,3,2,2)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSFRFF/ ALPFFQ,AKAPPA,GMUFFQ,SWEFF2
C
      IF (LPAR(2).EQ.0) THEN
        RA1=0D0
        RA3=0D0
        RB1=0D0
        RB3=0D0
        RETURN
      ENDIF
C
      LQF=MOD(LQ,2)
      IF (LQF.EQ.0) THEN
        LQF=3
        CQ=-1D0/3D0
        ELSE
        LQF=2
        CQ=2D0/3D0
      ENDIF
C
      T=-(SP-MEI2-MPRO2)*X*Y
      DG=1D0/T
      DZ=1D0/(T-MZ2)
C
C..COUPLING CONSTANTS INCLUDING RUNNING ALPHA AND WEAK FORM FACTORS
      CALL HSSETF(T)
      VEZ=VAFI1(1,1,2)
      AEZ=VAFI1(2,1,2)
C
C..STRUCTURE FUNCTIONS
      F1GG=FLIND1(1,LQF,1,1)
      F1GZ=FLIND1(1,LQF,1,2)
      F1ZZ=FLIND1(1,LQF,2,2)
      F3GG=-FLIND1(2,LQF,1,1)
      F3GZ=-FLIND1(2,LQF,1,2)
      F3ZZ=-FLIND1(2,LQF,2,2)
      EPEGG=FLIND1(1,1,1,1)
      EPEGZ=FLIND1(1,1,1,2)
      EPEZZ=FLIND1(1,1,2,2)
      EMEGG=-FLIND1(2,1,1,1)
      EMEGZ=-FLIND1(2,1,1,2)
      EMEZZ=-FLIND1(2,1,2,2)

C..SOFT BREMSSTRAHLUNG (INCLUDING QUARK SELF ENERGY)
      IF (LPAR(11).EQ.1) THEN
        BSQ=HSBRNC(X,Y,LL,LQ)
        RBSQ=DREAL(BSQ)*ALP2PI
        CALL HSSAB0(X,Y,POL,LQ,RA10,RA30,RB10,RB30)
C..QED VERTEX
        LPRK15=LPAR(15)
        LPAR(15)=0
        CVEG=HSFHFB(T,1,1,1,ME2)
        LPAR(15)=LPRK15
        DVERTX=DREAL(CVEG)
        RBSQ=RBSQ+2D0*DVERTX
        ELSE
        RBSQ=0D0
        RA10=0D0
        RA30=0D0
        RB10=0D0
        RB30=0D0
      ENDIF

C..BOX DIAGRAMS
      IF (LPAR(14).EQ.1) THEN
        S=SP*X
        U=-(T+S)
        BXGG = -LL*ALP1PI*DREAL(HSBCGA(T,S)-HSBCGA(T,U))
        BXGG5=  ALP1PI*DREAL(HSBCGA(T,S)+HSBCGA(T,U))
        BXGZ = -LL*2D0*ALP1PI*DLOG((MZ2-T)/(-T))*DLOG(-U/S)
     *         -LL*ALP1PI*DREAL(HSBXCV(T,S,CMZ2)-HSBXCV(T,U,CMZ2))
        BXGZ5=  ALP1PI*DREAL( HSBXCA(T,S,CMZ2) + HSBXCA(T,U,CMZ2) )
        F1BGG=F1GG*CQ
        F1BGZ=F1GZ*CQ
        F1BZZ=F1ZZ*CQ
        F3BGG=F3GG*CQ
        F3BGZ=F3GZ*CQ
        F3BZZ=F3ZZ*CQ
      ELSE
        BXGG=0D0
        BXGG5=0D0
        BXGZ=0D0
        BXGZ5=0D0
        F1BGG=0D0
        F1BGZ=0D0
        F1BZZ=0D0
        F3BGG=0D0
        F3BGZ=0D0
        F3BZZ=0D0
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        S=SP*X
        U=-(T+S)
        F1BZZ=F1ZZ*CQ
        F3BZZ=F3ZZ*CQ
        BXZZ = -LL*ALP1PI*DREAL(HSBXI0(T,U,CMZ2) - HSBXI0(T,S,CMZ2))
        BXZZ5=     ALP1PI*DREAL(HSBXI5(T,U,CMZ2) + HSBXI5(T,S,CMZ2))
        FGBWW = (VAFI1(1,LQF,1)+VAFI1(2,LQF,1))
     *         *(VAFI1(1,  1,1)+VAFI1(2,  1,1))/4D0/SW2/SW2
        FZBWW = (VAFI1(1,LQF,2)+VAFI1(2,LQF,2))
     *         *(VAFI1(1,  1,2)+VAFI1(2,  1,2))/4D0/SW2/SW2
        FBZZZ= FLIND1(1,LQF,2,2)*VAFI1(1,LQF,2)
     *        +FLIND1(2,LQF,2,2)*VAFI1(2,LQF,2)
        FBZZ5= FLIND1(1,LQF,2,2)*VAFI1(2,LQF,2)
     *        +FLIND1(2,LQF,2,2)*VAFI1(1,LQF,2)
        IF (LL.LT.0) THEN
          IF (MOD(LQ,2).EQ.1) THEN
            BXWW =  ALP1PI*DREAL(-HSBXI0(T,S,CMW2)+HSBXI5(T,S,CMW2))
            ELSE
            BXWW =  ALP1PI*DREAL( HSBXI0(T,U,CMW2)+HSBXI5(T,U,CMW2))
          ENDIF
        ELSE
          IF (MOD(LQ,2).EQ.1) THEN
            BXWW =  ALP1PI*DREAL(-HSBXI0(T,U,CMW2)+HSBXI5(T,U,CMW2))
            ELSE
            BXWW =  ALP1PI*DREAL( HSBXI0(T,S,CMW2)+HSBXI5(T,S,CMW2))
          ENDIF
        ENDIF
      ELSE
        BXZZ=0D0
        BXZZ5=0D0
        FGBWW=0D0
        FZBWW=0D0
        FBZZZ=0D0
        FBZZ5=0D0
        BXWW=0D0
      ENDIF
C
C..CONSTRUCTING THE CROSS SECTION
      IF (LPAR(3).LT.3) THEN
        A1=F1GG*DG*DG*EPEGG+2D0*F1GZ*DG*DZ*EPEGZ+F1ZZ*DZ*DZ*EPEZZ
     *     +RBSQ*RA10
        ELSE
        A1=(F1GG*DG*DG*EPEGG+2D0*F1GZ*DG*DZ*EPEGZ+F1ZZ*DZ*DZ*EPEZZ)
     *     *(1D0+RBSQ)
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        A1 = A1 +
     *     BXGG*EPEGG*DG*DG*F1BGG
     *    +(BXGG+BXGZ)*EPEGZ*DG*DZ*F1BGZ
     *    -(BXGG5+BXGZ5)*EMEGZ*DG*DZ*F3BGZ
     *    +BXGZ*EPEZZ*DZ*DZ*F1BZZ - BXGZ5*EMEZZ*DZ*DZ*F3BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        A1 = A1 +
     *     BXZZ*FLIND1(1,1,2,2)*DG*DG*F1BZZ
     *    -BXZZ5*DG*DG*FLIND1(2,1,2,2)*F3BZZ
     *    +BXZZ*(FLIND1(1,1,2,2)*VAFI1(1,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(2,1,2))*DG*DZ*FBZZZ
     *    -BXZZ5*(FLIND1(1,1,2,2)*VAFI1(2,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(1,1,2))*DG*DZ*FBZZ5
     *    +DG*(DG*FGBWW+DZ*FZBWW)*BXWW
      ENDIF
C
      IF (LPAR(3).LT.3) THEN
        B3=F3GG*DG*DG*EMEGG+2D0*F3GZ*DG*DZ*EMEGZ+F3ZZ*DZ*DZ*EMEZZ
     *     +RBSQ*RB30
        ELSE
        B3=(F3GG*DG*DG*EMEGG+2D0*F3GZ*DG*DZ*EMEGZ+F3ZZ*DZ*DZ*EMEZZ)
     *     *(1D0+RBSQ)
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        B3 = B3
     *    -EPEGG*DG*DG*F1BGG*BXGG5
     *    +(BXGG+BXGZ)*EMEGZ*DG*DZ*F3BGZ
     *    -(BXGG5+BXGZ5)*EPEGZ*DG*DZ*F1BGZ
     *    +BXGZ*EMEZZ*DZ*DZ*F3BZZ - BXGZ5*EPEZZ*DZ*DZ*F1BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        B3 = B3 +
     *     BXZZ*FLIND1(2,1,2,2)*DG*DG*F3BZZ
     *    -BXZZ5*DG*DG*FLIND1(1,1,2,2)*F1BZZ
     *    +BXZZ*(FLIND1(1,1,2,2)*VAFI1(2,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(1,1,2))*DG*DZ*FBZZ5
     *    -BXZZ5*(FLIND1(1,1,2,2)*VAFI1(1,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(2,1,2))*DG*DZ*FBZZZ
     *    +DG*(DG*FGBWW+DZ*FZBWW)*BXWW
      ENDIF
      RA1=A1
      RB3=B3
C
      IF (POL.EQ.0D0) THEN
        RB1=0D0
        RA3=0D0
        ELSE
        IF (LPAR(3).LT.3) THEN
          B1=F1GG*DG*DG*EMEGG+2D0*F1GZ*DG*DZ*EMEGZ+F1ZZ*DZ*DZ*EMEZZ
     *       +RBSQ*RB10
          ELSE
          B1=(F1GG*DG*DG*EMEGG+2D0*F1GZ*DG*DZ*EMEGZ+F1ZZ*DZ*DZ*EMEZZ)
     *       *(1D0+RBSQ)
        ENDIF
C
        IF (LPAR(14).EQ.1) THEN
          B1 = B1
     *      +(BXGG+BXGZ)*EMEGZ*DG*DZ*F1BGZ
     *      -(BXGG5+BXGZ5)*EPEGZ*DG*DZ*F3BGZ
     *      +BXGZ*EMEZZ*DZ*DZ*F1BZZ - BXGZ5*EPEZZ*DZ*DZ*F3BZZ
        ENDIF
        IF (LPAR(16).EQ.1) THEN
          B1 = B1 +
     *       BXZZ*EMEZZ*DG*DG*F1BZZ - BXZZ5*DG*DG*EPEZZ*F3BZZ
     *      -BXZZ*(EPEZZ*AEZ-EMEZZ*VEZ)*DG*DZ*FBZZZ
     *      +BXZZ5*(EPEZZ*VEZ-EMEZZ*AEZ)*DG*DZ*FBZZ5
     *      -DG*(DG*FGBWW+DZ*FZBWW)*BXWW
        ENDIF
C
        IF (LPAR(3).LT.3) THEN
          A3=F3GG*DG*DG*EPEGG+2D0*F3GZ*DG*DZ*EPEGZ+F3ZZ*DZ*DZ*EPEZZ
     *       +RBSQ*RA30
          ELSE
          A3=(F3GG*DG*DG*EPEGG+2D0*F3GZ*DG*DZ*EPEGZ+F3ZZ*DZ*DZ*EPEZZ)
     *       *(1D0+RBSQ)
        ENDIF
C
        IF (LPAR(14).EQ.1) THEN
          A3 = A3 +
     *       (BXGG+BXGZ)*EPEGZ*DG*DZ*F3BGZ
     *      -(BXGG5+BXGZ5)*EMEGZ*DG*DZ*F1BGZ
     *      +BXGZ*EPEZZ*DZ*DZ*F3BZZ - BXGZ5*EMEZZ*DZ*DZ*F1BZZ
        ENDIF
        IF (LPAR(16).EQ.1) THEN
          A3 = A3 +
     *       BXZZ*EPEZZ*DG*DG*F3BZZ - BXZZ5*DG*DG*EMEZZ*F1BZZ
     *      +BXZZ*(EPEZZ*VEZ-EMEZZ*AEZ)*DG*DZ*FBZZ5
     *      -BXZZ5*(EPEZZ*AEZ-EMEZZ*VEZ)*DG*DZ*FBZZZ
     *      -DG*(DG*FGBWW+DZ*FZBWW)*BXWW
        ENDIF
        RB1=B1
        RA3=A3
      ENDIF
      RETURN
      END
