C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSNCG2(XARG)
C***
C   LAST CHANGE 24/07/90 BY HJM
C***
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DOUBLE PRECISION HSNCG2
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      DIMENSION XARG(2)
C
      DX=XMAX-XMIN
      X=XMIN+XARG(1)*DX
      Z=XARG(2)
      GL=HSGLOW(X)
      GU=HSGUPP(X)
      DG=DMAX1(GU-GL,0D0)
      G=GL+Z*DG
      Q2=-1D0/G
      IF(IPRINT.GT.20) WRITE(LUNTES,'(A,4D15.6)')
     &                      ' HSNCG2: X, Z, G, Q2',X,Z,G,Q2
      HSNCG2=Q2**2*HSNC22(X,Q2)*DG*DX
      RETURN
      END
