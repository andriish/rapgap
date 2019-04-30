      
      INTEGER FUNCTION ICON_PDG_SIB(ID)
C************************************************************************
C
C   convert PDG particle codes to SIBYLL particle codes
C
C                                         (R.E. 09/97)
C
C************************************************************************
      SAVE

      DIMENSION ITABLE(49)
      DATA ITABLE /
*  gam
     &                  22,
*  e+
     &                 -11,
*  e-
     &                  11,
*  mu+
     &                 -13,
*  mu-
     &                  13,
*  pi0
     &                 111,
*  pi+
     &                 211,
*  pi-
     &                -211,
*  k+
     &                 321,
*  k-
     &                -321,
*  k0l
     &                 130,
*  k0s
     &                 310,
*  p
     &                2212,
*  n
     &                2112,
*  nue
     &                  12,
*  nueb
     &                 -12,
*  num
     &                  14,
*  numb
     &                 -14,
*  pbar
     &           -99999999,
*  nbar
     &           -99999999,
*  k0
     &                 311,
*  k0b
     &                -311,
*  eta
     &                 221,
*  etap
CHS     &                 311,
CHS
CHS                       should read 331
CHS changed:
     &                 331,
*  rho+
     &                 213,
*  rho-
     &                -213,
*  rho0
     &                 113,
*  k*+
     &                 323,
*  k*-
     &                -323,
*  k*0
     &                 313,
*  k*0b
     &                -313,
*  omeg
     &                 223,
*  phi
     &                 333,
*  SIG+
     &                3222,
*  SIG0
     &                3212,
*  SIG-
     &                3112,
*  XI0
     &                3322,
*  XI-
     &                3312,
*  LAM
     &                3122,
*  DELT++
     &                2224,
*  DELT+
     &                2214,
*  DELT0
     &                2114,
*  DELT-
     &                1114,
*  SIG*+
     &                3224,
*  SIG*0
     &                3214,
*  SIG*-
     &                3114,
*  XI*0
     &                3324,
*  XI*-
     &                3314,
*  OME*-
     &                3334 /



      IDPDG = ID

 100  CONTINUE
      IDA = ABS(ID)

      IF(IDA.GT.1000) THEN
        IS = IDA
        IC = SIGN(1,IDPDG)
      ELSE
        IS = IDPDG
        IC = 1
      ENDIF

      DO I=1,49
        IF(IS.EQ.ITABLE(I)) THEN
          ICON_PDG_SIB = I*IC
          RETURN
        ENDIF
      ENDDO

      IF(IDPDG.EQ.80000) THEN
        ICON_PDG_SIB = 13
      ELSE  
*       IF(IDA.EQ.411) THEN
*         IDPDG = SIGN(311,IDPDG)
*       ELSE IF(IDA.EQ.421) THEN
*         IDPDG = SIGN(321,IDPDG)
*       ELSE IF(IDA.EQ.431) THEN
*         IDPDG = SIGN(331,IDPDG)
*       ELSE IF(IDA.EQ.4122) THEN
*         IDPDG = SIGN(3122,IDPDG)
*       ELSE IF(IDA.EQ.4222) THEN
*         IDPDG = SIGN(3222,IDPDG)
*       ELSE IF(IDA.EQ.4222) THEN
*         IDPDG = SIGN(3222,IDPDG)
*       ELSE
          print *,'ICON_PDG_DTU: no particle found for ',IDPDG
          ICON_PDG_SIB = 0
          RETURN
*       ENDIF
*       GOTO 100
      ENDIF

      END
