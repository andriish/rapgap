C
      SUBROUTINE DJGVAR(ICHNN,X,Y,Q2)
C---
C   TRANSFER VARIABLES FROM HERACLES TO LEPTO OR OTHER ROUTINE
C   FOR FRAGMENTATION
C---
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
#include "hepevt.inc"
      COMMON /HSRESC/ SSCH,Q2SCH,W2SCH,XSCH,YSCH
     &               ,SLIX,SLIY,SLIZ,SLIE,SLIM
     &               ,SLFX,SLFY,SLFZ,SLFE,SLFM

C---DECLARATIONS FOR LEPTO
C---->
C   TRANSER OF KINEMATICAL VARIABLES
C---->
      B=X
      RETURN
      END
