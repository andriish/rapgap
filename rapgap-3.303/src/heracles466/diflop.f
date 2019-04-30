
C*********************************************************************
C...Translate input from FLONG input flag

      SUBROUTINE DIFLOP
      COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U

      PARL(11)=PARL11
      PARL(19)=PARL19
      RETURN
      END
