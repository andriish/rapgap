      SUBROUTINE DIFFR4(X,WDIF)
      IMPLICIT NONE
#include "rgparas.inc"
      Double Precision GEV2NB,WPART,WDIF,WT1
	Double Precision X(20)
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      WT1=0.D0
      CALL PARTDF(X,WPART)
c      write(6,*) 'diffr4: WPART ',WPART
      IF(WPART.NE.0.) THEN
         IF(IPRO.EQ.100) CALL ELERHO(WT1)
      ENDIF
      WDIF=WPART*WT1*GEV2NB
c      write(6,*) 'diffr4: WDIF ',WDIF
c      call lulist(1)
      RETURN
      END


