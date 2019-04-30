C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCWQR(FS,CM2,MQ2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSCWQR,HSSPEN,HSCLN,CSX,CM2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      CSX=DCMPLX(FS,1D-6)
      HSCWQR=-1D0/FS*(HSCLN(-CSX/CM2)*HSCLN(-CSX/CM2)/2D0
     .                +HSSPEN((CSX-CM2)/CSX)-PI*PI/3D0
     .                +DCMPLX(0D0,1D0)*PI*HSCLN((CSX-CM2)/CSX))
      RETURN
      END
