C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     INTEGRATION SUBROUTINE
C     SUBSTITUTE FOR THE NAGLIB ROUTINE D01FCF
C     Note: The function to be integrated is identified by 
C     INTEGER IFUN = 1: Neutral Current
C                  = 2: Charged Current
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DX1FCF(NDIMEN,BLOW,BUP,MINPTS,MAXPTS,IFUN,ACCREQ,
     &              ACCFIN,LENWRK,WRKSTR,RESULT,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BLOW(2),BUP(2),BLOWC(2),BUPC(2),WRKSTR(LENWRK)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSD01L/ BLOWC,BUPC,ARG1,ACC,INCCC
      EXTERNAL DFNCII
      DO 1 I=1,2
      BLOWC(I)=BLOW(I)
    1 BUPC(I)=BUP(I)
      ACC=ACCREQ
      INCCC=IFUN
      RESULT=GAUSK1(DFNCII,BLOWC(2),BUPC(2),ACCREQ)
      ACCFIN=ACCREQ
      IFAIL=0
      RETURN
      END
