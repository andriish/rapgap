C######################################################################
C
C   Various routines to give structure function parametrizations.
C   All but LNSTRF can be used separately without initialization.
C
C ********************************************************************

      SUBROUTINE LNSTRF(X,Q2,XPQ)

C...Structure function per nucleon for a proton/neutron mixture
C...according to defined nucleus.

      DIMENSION XPQ(-6:6)

      CALL RGSTFU(2212,X,Q2,XPQ)

      RETURN
      END
