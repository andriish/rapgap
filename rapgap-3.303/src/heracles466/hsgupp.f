C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSGUPP(X)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSPCUT/ PTMIN,PTXM0
C
      GSP=SP-MEI2-MPRO2
      YMAXX=X*(1D0-4D0*MEI2*MPRO2/GSP/GSP)
     *        /(X*(1D0+X*MPRO2/GSP)+MEI2/GSP)
      IF(ICUT.EQ.1) THEN
C                                   CUT IN EXTERNALLY DEFINED Q**2(MIN)
        GMAX=-1D0/(X*YMAXX*GSP)
      ELSEIF(ICUT.EQ.2) THEN
C                                   CUT IN W / MAXIMUM Q**2
        GMAX=-1D0/(X*YMAXX*GSP)
      ELSEIF(ICUT.EQ.3) THEN
C                                   CUT IN Y
        Q2MAX1=X*MIN(YMAX,YMAXX)*GSP
C                                   CUT ON ELECTRON TRANSVERSE MOMENTUM
C                                   (MASSES NEGLECTED)
        QP2MAX=X*GSP/2D0*(1D0+DSQRT(1D0-PTXM0/X))
        GMAX=-1D0/MIN(Q2MAX1,QP2MAX,Q2MAX)
      ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSGUPP'
        STOP
      ENDIF
C
      HSGUPP=GMAX
      IF(HSGUPP.LT.GMIN) THEN
        HSGUPP=GMIN
      ENDIF
      RETURN
      END
