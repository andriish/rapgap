C->
      BLOCK DATA PARAM_INI
C....This block data contains default values
C.   of the parameters used in fragmentation
C................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      SAVE
      COMMON /S_CZDIS/ FA, FB0
      COMMON /S_CZDISs/ FAs1, fAs2
      COMMON /S_CZLEAD/ CLEAD, FLEAD
      COMMON /S_CPSPL/ CCHIK(3,6:14)
      COMMON /S_CQDIS/ PPT0 (33),ptflag
      COMMON /S_CDIF0/ FFD, FBD, FDD
      COMMON /S_CFLAFR/ PAR(8)
C...Longitudinal Fragmentation function
      DATA FA /0.5/, FB0 /0.8/
C...Longitudinal Fragmentation function for leading baryons
       DATA CLEAD  /0.0/, FLEAD  /0.6/
c      strange fragmentation
      data FAs1 /3./, fAs2 /3./
c      data FAs1 /0./, fAs2 /0./
C...pT of sea partons
      DATA PTFLAG /1./
      DATA PPT0 /0.30,0.30,0.450,30*0.60/
C...Splitting parameters
      DATA CCHIK /21*2.,6*3./
C...Parameters of flavor formation
      DATA PAR /0.04,0.25,0.25,0.14,0.3,0.3,0.15,0./
      END
