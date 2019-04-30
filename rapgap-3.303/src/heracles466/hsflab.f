C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   PARAMETERS OF THE HERA-LAB-SYSTEM REAL*8
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFLAB(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      EH = EELE
      PH = PELE
      EQH = XS*EPRO
      PQH = XS*PPRO
      ESH = -(PPRO*(TP-MEI2-MEF2)+PH*(UP-MPRO2-MEF2))
     &              /2D0/(EH*PPRO+PH*EPRO)
      PSH = DSQRT((ESH+MEF)*(ESH-MEF))
      COSEH = (TP-2D0*MEF2+2D0*EH*ESH)/2D0/PH/PSH
      SINH2=1D0-COSEH*COSEH
      IF(SINH2.LE.0D0) THEN
        SINEH=0D0
        ELSE
        SINEH=DSQRT(SINH2)
      ENDIF
      END
