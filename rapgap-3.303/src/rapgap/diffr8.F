      SUBROUTINE DIFFR8(X,WDIF)
      IMPLICIT None
	Double Precision X,WDIF,WPART,WT1,GEV2NB
#include "rgparas.inc"
#include "rghs45.inc"
      DIMENSION X(20)
	Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      DATA GEV2NB/.3893857D+6/
      WPART=0.D0
      WDIF=0.D0
      WT1=0.D0
      IF(IHERAC.EQ.0) THEN
         CALL PARTDF(X,WPART)
      ELSE
         write(6,*) ' resloved photon with HERACLES not implemented '
      ENDIF
c      write(6,*) ' DIFFR8 WPART before eleres ',WPART,WT1
      IF(WPART.GT.0.) THEN
         CALL ELERES(WT1)
ccc	      wt1=1
      ENDIF
c      write(6,*) ' DIFFR7 after eleres WPART,WT1 ',WPART,WT1
      WDIF=WPART*WT1*GEV2NB
      IF(WDIF.LE.0.D0.AND.IDEBUG.EQ.1) THEN
         write(6,*) ' DIFFR8 WPART,WT1 ',WPART,WT1
      ENDIF
c      IF(WDIF.LE.0.D0) THEN
c         write(6,*) ' DIFFR8 WPART,WT1 ',WPART,WT1
c      ENDIF
      RETURN
      END
