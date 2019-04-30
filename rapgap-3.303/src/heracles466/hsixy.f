C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSIXY(X,Y)
      COMPLEX*16 HSIXY,X,Y,HSSPEN,HSCLN
      HSIXY=HSSPEN((1D0-Y)/(X-Y))-HSSPEN(Y/(Y-X))
     *    + HSCLN((1D0-X)/(Y-X))*(HSCLN(1D0-Y)-HSCLN(X-Y))
     *    - HSCLN(X/(X-Y))*(HSCLN(-Y)-HSCLN(X-Y))
     *    + HSCLN((X-1D0)/X)*HSCLN(X-Y)
      END
