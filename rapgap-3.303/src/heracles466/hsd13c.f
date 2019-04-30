C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSD13C(FS,FT,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX*16 HSD13C,CM2,HSSPEN,HSCLN
      FU=-(FS+FT)
      HSD13C=1D0/FU*(HSSPEN(-FS/CM2)-HSSPEN(FT/CM2)
     *               +HSCLN(-FS/CM2)*HSCLN((CM2+FS)/(CM2-FT)))
      END
