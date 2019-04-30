      SUBROUTINE PARTHS(X,WPART)
      IMPLICIT NONE
#include "rgdisdif.inc"
#include "rglujets.inc"
	Double PRecision X,WPART
      DIMENSION X(20)
      IF(IDIR.EQ.0) THEN
         CALL PARTDFHS(X,WPART)
      ELSE
         CALL PARTDIHS(X,WPART)
      ENDIF
      RETURN
      END
