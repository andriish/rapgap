C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     INTEGRATION SUBROUTINE
C     SUBSTITUTE FOR THE NAGLIB ROUTINE D01AJF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE D01AJF(F,BLOW,BUP,ACCABS,ACCREL,RESULT,ACCFIN,
     *                W,LW,IW,LIW,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IW(LIW),W(LW)
C      EXTERNAL HSELG1
      EXTERNAL F
      BLOWC=BLOW
      BUPC=BUP
      ACCREQ=ACCREL
C      RESULT=GAUSK1(HSELG1,BLOWC,BUPC,ACCREQ)
      RESULT=GAUSK1(F,BLOWC,BUPC,ACCREQ)
      ACCFIN=ACCREL*RESULT
      IFAIL=0
      RETURN
      END
