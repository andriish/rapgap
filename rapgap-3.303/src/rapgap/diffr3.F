      SUBROUTINE DIFFR3(X,WDIF)
      IMPLICIT NONE
#include "rghs45.inc"
      Double Precision GEV2NB,WPART,WDIF
	Double Precision X(20)
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      IF(IHERAC.EQ.0) THEN
         CALL PARTDF(X,WPART)
      ELSE
         CALL PARTDFHS(X,WPART)
      ENDIF
      IF(WPART.LE.0.D0) WPART = 0.D0
c matrix element in partdf included
      WDIF=WPART*GEV2NB
      RETURN
      END
