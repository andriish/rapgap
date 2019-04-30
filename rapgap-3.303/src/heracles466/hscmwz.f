C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCMWZ(GS,CM1,CM2)
      DOUBLE PRECISION GS
      COMPLEX*16 HSCMWZ,CM1,CM2,X1,X2,C,HSCLN
      C=CDSQRT((GS-CM1-CM2)*(GS-CM1-CM2)-4D0*CM1*CM2)
      X1=(GS-CM2+CM1+C)/2D0/GS
      X2=(GS-CM2+CM1-C)/2D0/GS
      HSCMWZ=-1D0/GS*HSCLN(X1/(X1-1D0))*HSCLN(X2/(X2-1D0))
      END
