      SUBROUTINE DIR3(X,WDIF)
      IMPLICIT NONE
#include "rghs45.inc"
      Double Precision GEV2NB,WPART,WDIF
      Double Precision X(20)
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      IF(IHERAC.EQ.0) THEN
         CALL PARTDI(X,WPART)
      ELSE
         CALL PARTDIHS(X,WPART)
      ENDIF
c matrix element in partdi included
      IF(WPART.LE.0.D0) WPART=0.D0
      WDIF=WPART*GEV2NB
c      write(6,*) ' dir3 ',wpart,wdif
      RETURN
      END
