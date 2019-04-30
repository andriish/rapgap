      SUBROUTINE DIR1(X,WDIF)
      IMPLICIT NONE
#include "rgparas.inc"
#include "rghs45.inc"
#include "rgnlo.inc"
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
	   IF(IBS.EQ.0) THEN
            IF(IPRO.EQ.10) CALL ELEQQL(WT1)
            IF(IPRO.EQ.13) CALL ELEQQF(WT1)
            IF(IPRO.EQ.15.OR.IPRO.EQ.16) CALL ELEQCDC(WT1)
	      ELSEIF(IBS.EQ.1) THEN
            IF(IPRO.EQ.10) CALL ELQQLNL(WT1)
            IF(IPRO.EQ.13) CALL ELQQFNL(WT1)
            IF(IPRO.EQ.15) CALL ELQCDCNL(WT1)
		ELSE
		write(6,*) ' IBS not implemented ',IBS
	   ENDIF
         WDIF=WPART*WT1*GEV2NB
      ENDIF
c           write(6,*) 'dir1 ',wpart,wt1,IHERAC
      RETURN
      END
