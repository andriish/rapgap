

      INTEGER FUNCTION ICON_SIB_PDG(ISIB)
C************************************************************************
C
C   convert SIBYLL particle codes to PDG particle codes
C
C                                         (H.S. 07/99)
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

      ISIBA=ABS(ISIB)
      IF (ISIBA.GE.50) THEN
        print *,'ICON_SIB_PDG: wrong particle code: ISIB= ',ISIB
        ICON_SIB_PDG=0
        RETURN
      ENDIF
      ISGN=SIGN(1,ISIB)
      ICON_SIB_PDG=ISGN*ITABLE(ISIBA)
      RETURN
      END
