C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCMW(FS,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSCMW,HSSPEN,CM2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      HSCMW=-1D0/FS*(HSSPEN((FS+CM2)/CM2)-PI*PI/6D0)
      END
