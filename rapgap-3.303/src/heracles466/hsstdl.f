
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   Interface to the parametrization of DL for F2
C   A.Donnachie and P.Landshoff, DAMTP 93-23 and M/C-TH 93/11
C   Published in Z.Phys.C61 (1994) (hep-ph/9305319) 

      SUBROUTINE HSSTDL(X,Q2,ZF1,ZF2)
      DOUBLE PRECISION X,Q2,ZF1,ZF2,HSF2DL
      ZF2=HSF2DL(Q2,X)
      ZF1=ZF2/2D0/X
      RETURN
      END
