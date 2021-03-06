C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSELG2(XARG)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      DIMENSION XARG(1)
C
      X=1D0
      Z=XARG(1)
      GSP=SP-MEI2-MPRO2
      Q2MNY=YMIN*GSP
      GL=-1D0/MAX(Q2MIN,Q2MNY)
      YMAXX=(1D0-4D0*MEI2*MPRO2/GSP/GSP)/(1D0+X*MPRO2/GSP+MEI2/GSP)
      GU=-1D0/(MIN(YMAX,YMAXX)*GSP)
      DG=GU-GL
      G=GL+Z*DG
      Q2=-1D0/G
      IF(IPRINT.GE.20) WRITE(LUNTES,'(A,3D15.6)')
     &                      ' HSELG2: Z, G, Q2',Z,G,Q2
      HSELG2=Q2**2*HSEL22(Q2)*DG
      RETURN
      END
