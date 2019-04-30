C **********************************************************************

      FUNCTION FLQINT(Z)

C...Quark contribution integrand to QCD longitudinal structure function.
      DOUBLE PRECISION PARL
      COMMON /RAPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
#include "rgluco.inc"
      DATA PI/3.14159/
      NTOT=NTOT+1
      IF(IHF.EQ.0) THEN
         CALL LNSTRF(Z,Q2,XPQ)
         FLQINT=0.
         DO 10 I=-LST(12),LST(12)
            IF(I.EQ.0) GOTO 10
            FLQINT=FLQINT+QC(IABS(I))**2*XPQ(I)
   10    CONTINUE
         FLQINT=4./3.*PARL(25)/PI*(X/Z)**2*FLQINT/Z
      ENDIF
      NPASS=NPASS+1


      RETURN
      END
