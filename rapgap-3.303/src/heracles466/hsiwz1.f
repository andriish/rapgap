C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       FUNCTION HSIWZ1(GS,GT,CM1,CM2)
       COMPLEX*16 HSIWZ1,HSD0,HSCMWZ,CM1,CM2
       DOUBLE PRECISION GS,GT
       HSIWZ1=GT/2D0*HSD0(GS,GT,CM1,CM2)-HSCMWZ(GS,CM1,CM2)
       END
