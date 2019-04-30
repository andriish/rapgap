C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND FOR NONRADIATIVE CONTRIBUTION IN CHARGED CURRENT
C   (ARGUMENTS AS REQUIRED BY INTEGRATION ROUTINE IN INIT22;
C    HSCCG1 AND HSCCG2 DIFFER ONLY IN THE LIST OF ARGUMENTS)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCG1(NDIMEN,ARGUM)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSINTL/ XL,XU
      DIMENSION ARGUM(NDIMEN)
C
      DX=XU-XL
      X=XL+ARGUM(1)*DX
      GS=SP-MEI2-MPRO2
      YMAXX=X*(1D0-4D0*MEI2*MPRO2/GS/GS)/(X*(1D0+X*MPRO2/GS)+MEI2/GS)
      Q2L=Q2MIN
      Q2U=X*GS
      IF(ICUT.EQ.2) THEN
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        Q2L=MAX(Q2MIN,QQ2MIN)
        Q2U=X*GS
      ELSEIF(ICUT.GE.3) THEN
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        QQQ2MN=X*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+X*CTHCON)
        QT2MIN=X*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=X*GS/2D0*(1D0-DSQRT(1D0-PTXM0/X))
        QP2MAX=X*GS/2D0*(1D0+DSQRT(1D0-PTXM0/X))
        YP2MAX=QP2MAX/GS/X
        Q2L=MAX(Q2MIN,QQ2MIN,QQQ2MN,QT2MIN,QP2MIN)
        Q2U=X*MIN(YMAX,YMAXX,YP2MAX)*GS
        Q2U=MIN(Q2U,Q2MAX)
      ENDIF
      DQ2=DMAX1(Q2U-Q2L,0D0)
      Q2=Q2L+DQ2*ARGUM(2)
      IF(IPRINT.GT.20) WRITE(LUNTES,'(A,2D15.6)')
     &                      ' HSCCG1: X, Q2',X,Q2
      HSCCG1=HSCC22(X,Q2)*DX*DQ2
      RETURN
      END
