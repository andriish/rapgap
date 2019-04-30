
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION HSHADQ(S)
C  HADRONIC IRREDUCIBLE QQ SELF-ENERGY: TRANSVERSE
* THIS IS A SLIGHTLY MODIFIED VERSION OF BURKHARDTS CODE. DB. TR.,01/91
C     parametrize the real part of the photon self energy function
C     by  a + b ln(1+C*:S:) , as in my 1981 TASSO note but using
C     updated values, extended using RQCD up to 100 TeV
C     for details see:
C     H.Burkhardt, F.Jegerlehner, G.Penso and C.Verzegnassi
C     in CERN Yellow Report on "Polarization at LEP" 1988
C     H.BURKHARDT, CERN/ALEPH, AUGUST 1988
C     negative values mean t - channel (spacelike)
C     positive values mean s - channel (timelike )
C     in the space like values around 1 GeV are typical for luminosity
C     the values at 92 GeV ( Z mass ) give the light quark contribution
C     to delta r
C     take care of the sign of REPI when using this in different
C     programs
C     Here REPI was chosen to
C     be positive (so that it corresponds directly to delta alpha)
C     often its assumed to be negative.
C
C     the imaginary part is proportional to R (had / mu cross section)
C     and is therefore 0 below threshold ( in all the spacelike region)
C     note also that alpha_s usually has been derived from the measured
C     values of R.
C     Changing alpha_s to values incompatible with current data
C     would imply to be also inconsistent with RE,IM PI
C     defined here
C
C     H.BURKHARDT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DATA A1,B1,C1/   0.0   ,  0.00835,  1.0  /
      DATA A2,B2,C2/   0.0   ,  0.00238,  3.927 /
      DATA A3,B3,C3/ 0.00165 ,  0.00300,  1.0  /
      DATA A4,B4,C4/ 0.00221 ,  0.00293,  1.0  /
C
      DATA PI/3.141592653589793D0/,ALFAIN/137.0359895D0/,INIT/0/
C
      IF(INIT.EQ.0) THEN
        INIT=1
        ALFA=1./ALFAIN
        ALFAPI=1./PI/ALFAIN
      ENDIF
      T=ABS(S)
      IF(T.LT.0.3**2) THEN
        REPIAA=A1+B1*LOG(1.+C1*T)
      ELSEIF(T.LT.3.**2) THEN
        REPIAA=A2+B2*LOG(1.+C2*T)
      ELSEIF(T.LT.100.**2) THEN
        REPIAA=A3+B3*LOG(1.+C3*T)
      ELSE
        REPIAA=A4+B4*LOG(1.+C4*T)
      ENDIF
C     as imaginary part take -i alfa/3 Rexp
      HSHADQ=REPIAA
      END
