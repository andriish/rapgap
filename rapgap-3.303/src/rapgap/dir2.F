      SUBROUTINE DIR2(X,WDIF)
      IMPLICIT NONE
#include "rgparas.inc"
#include "rghs45.inc"
      Double Precision GEV2NB,WPART,WDIF,WT1
      Double Precision X(20)
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      WT1=0.D0
      IF(IHERAC.EQ.0) THEN
         CALL PARTDI(X,WPART)
      ELSE
         CALL PARTDIHS(X,WPART)
      ENDIF

      IF(WPART.GT.0.) THEN
         IF(IPRO.EQ.11) CALL ELEQQB(WT1)
         IF(IPRO.EQ.14) CALL ELEQQF(WT1)
      ENDIF
      WDIF=WPART*WT1*GEV2NB
c          write(6,*) 'dir2 ',wpart,wt1,iherac
      RETURN
      END
