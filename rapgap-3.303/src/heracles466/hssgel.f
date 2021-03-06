C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   CROSS SECTION FOR ELASTIC LEPTON PROTON SCATTERING:
CCCCCCC 01.12.95    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC H. SPIESBERGER CCC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSSGEL(Q2,LL,POL,LQ)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO

      SHAT=SP-MEI2-MPRO2
      X=1D0
      Y=Q2/SHAT
      T=-Q2
      SPN=SXNORM*SP
      CALL HSFIVC(X,Y)
      CALL HSDELX(X,Y)
      CALL HSFIEL(-T,F1EL,F2EL)
      DO 10 I=1,12
   10 CQP(I)=0D0
      HSSGNC=0D0
      R1=Y*Y*(1D0+2D0*MEI2/T)
      R2=(1D0-Y+MPRO2*T/SHAT/SHAT)/X
      CQP(12)=F1EL*R1+F2EL*R2
      HSSGEL=8D0*SPN*CQP(12)/T/T
      END
