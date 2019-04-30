      SUBROUTINE LPRIKT(S,PT,PHI)
	IMPLICIT NONE
      Double Precision draprn
	Real S,PT,PHI
      EXTERNAL draprn
C...Size (PT) and azimuthal angle (PHI) of primordial kt according
C...to a Gaussian distribution.

      PT=S*SQRT(-ALOG(Real(draprn())))
      PHI=6.2832*Real(draprn())
      RETURN
      END
