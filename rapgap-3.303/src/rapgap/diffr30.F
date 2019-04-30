      SUBROUTINE DIFFR30(X,WDIF)
      IMPLICIT NONE
#include "rgparas.inc"
      Double Precision WPART,WDIF
	Double Precision X(20)
      WPART=0.D0
      WDIF=0.D0
c      write(6,*) 'diffr30 WPART before rgsatrap'
      CALL RGSATRAP(X,WPART)
c      write(6,*) 'diffr30 WPART after rgsatrap'
c      WDIF=WPART*GEV2NB
      WDIF=WPART
c      write(6,*) 'diffr30: WDIF ',WDIF
c      call lulist(1)
      RETURN
      END


