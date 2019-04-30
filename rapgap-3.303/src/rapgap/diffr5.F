      SUBROUTINE DIFFR5(X,WDIF)
      IMPLICIT None
	Double Precision X,WDIF,WPART,WT1,GEV2NB
#include "rgparas.inc"
#include "rghs45.inc"
      DIMENSION X(20)
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      WT1=0.D0
      IF(IHERAC.EQ.0) THEN
         CALL PARTDH(X,WPART)
      ELSE
         CALL PARTDHHS(X,WPART)
      ENDIF
      IF(WPART.GT.0.) THEN
         CALL ELEQQG(WT1)
      ENDIF
      if(wt1.eq.0d0.and.wpart.ne.0.d0) then
c      write(6,*) 'diffr5 wpart,wt1',wpart,wt1
      endif
c      write(6,*) 'diffr5 wpart,wt1',wpart,wt1
      WDIF=WPART*WT1*GEV2NB
      RETURN
      END
