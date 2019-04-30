
C*********************************************************************
C...translate input from ALFAS input

      SUBROUTINE DIALFS
      COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U

      MSTU(111)=MST111
      MSTU(115)=MST115
      PARU(111)=PAR111
      PARU(112)=PAR112
      RETURN
      END
