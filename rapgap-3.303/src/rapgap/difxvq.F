      SUBROUTINE DIFxvq(X,WDIF)
      IMPLICIT NONE
#include "rghs45.inc"
      Double Precision GEV2NB,WPART,WDIF,WT1
      Double Precision X(20)
	INTEGER IDEBUG
      COMMON/ INTERN/IDEBUG
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      WT1=0.D0
c	IDEBUG = 1
      IF(IHERAC.EQ.0) THEN
         CALL PARTDF(X,WPART)
      ELSE
         write(6,*) ' W/Z with HERACLES not implemented '
      ENDIF

      IF(WPART.GT.0.) THEN
         CALL ELEVQG(WT1)
      ENDIF
cc      write(6,*) wpart,wt1
      WDIF=WPART*WT1*GEV2NB
      IF(WDIF.LE.0.D0.AND.IDEBUG.EQ.1) THEN
         write(6,*) ' DIFxvq WPART,WT1 ',WPART,WT1
      ENDIF
      RETURN
      END
