C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSGLOW(X)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
C
      GSP=SP-MEI2-MPRO2
      IF(ICUT.EQ.1) THEN
C                                   CUT IN EXTERNALLY DEFINED Q**2(MIN)
        HSGLOW=-1D0/Q2MIN
      ELSEIF(ICUT.EQ.2) THEN
C                                   CUT IN W / MAXIMUM Q**2
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        HSGLOW=-1D0/MAX(Q2MIN,QQ2MIN)
      ELSEIF(ICUT.EQ.3) THEN
C                                   CUT IN Y AND W
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        QQQ2MN=X*YMIN*GSP
C                                   CUT ON ELECTRON SCATTERING ANGLE
C                                   (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+X*CTHCON)
        QT2MIN=X*YMINTH*GSP
C                                   CUT ON ELECTRON TRANSVERSE MOMENTUM
C                                   (MASSES NEGLECTED)
        QP2MIN=X*GSP/2D0*(1D0-DSQRT(1D0-PTXM0/X))
        HSGLOW=-1D0/MAX(Q2MIN,QQ2MIN,QQQ2MN,QT2MIN,QP2MIN)
      ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSGLOW'
        STOP
      ENDIF
      GMIN=HSGLOW
      RETURN
      END
