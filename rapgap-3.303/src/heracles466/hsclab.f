C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   PARAMETERS OF THE HERA-LAB-SYSTEM REAL*8 (FOR CC)
C   FROM CCLABEA (epcctot.f, 28/08/98)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSCLAB(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
C
      EH=EELE
      PH=PELE
      S=SP*XS
      CHI=(S-MEI2-MQI2)/2D0
      QH1=-(CHI*PELE+EELE*DSQRT(CHI*CHI-MEI2*MQI2))
      PQH=-(CHI*CHI-EELE*EELE*MQI2)/QH1
      EQH=DSQRT(PQH*PQH+MQI2)
      GS=SP-MEI2-MPRO2
      GU=-(1D0-Y)*GS
      TP=-X*Y*GS
      UP=GU+MPRO2
      ESH=-(PPRO*(TP-MEI2)+PH*(UP-MPRO2))/2D0/(EH*PPRO+PH*EPRO)
      PSH=ESH
      COSEH=(TP+2D0*EH*ESH)/2D0/PH/PSH
      SINH2=1D0-COSEH*COSEH
      IF(SINH2.LE.0D0) THEN
        SINEH=0D0
        ELSE
        SINEH=DSQRT(SINH2)
      ENDIF
      END
