C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   MINIMAL PHOTON ENERGY AS FUNCTION OF XS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSOMIN(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      CALL HSFLAB(X,Y,XS)
      STU = (XS-X)*Y*GS
      STUM = (XS-X)*Y*GS + XS*XS*MPRO2 + MEI2 + MEF2
      HO1 = STU/2D0/STUM*( EH + XS*EPRO - ESH
     *                     + DSQRT(PSH*PSH*SINEH*SINEH
     *                             + (PH - XS*PPRO - PSH*COSEH)
     *                              *(PH - XS*PPRO - PSH*COSEH)))
      HO2 = STU*STU/4D0/STUM/HO1
      HSOMIN = DMIN1(HO1,HO2)
      END
