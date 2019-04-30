      FUNCTION ALPHAS(RQ)
      IMPLICIT NONE
      Integer NF
      REAL ALAM
      COMMON/RGALPHA/NF,ALAM
#include "rgpara.inc"
#include "rgpypars.inc"
#include "rgdisdif.inc"
#include "rgludat1.inc"
#include "rgdiffr.inc"
#include "rgnlo.inc"
	Double Precision ALPHAS,RQ,RQ1,ALPHAS2
C PDFLIB COMMON
      Double Precision QCDL4,QCDL5
      COMMON/W50512/ QCDL4,QCDL5
C
      Double Precision SCALE
      RQ1 = RQ
      IF(IRUNA.EQ.0) THEN
         ALPHAS=ALPHS
      ELSE

         SCALE = RQ1*RQ1
	   
         ALPHAS =PYALPS(SCALE)
	   ALAM=PARU(112)
c	   write(6,*) ' ALPAHS old :', ALPHAS,scale,alam
	   IF(IBS.EQ.1) Then 
	     MSTU(111) = 2        
	     ALPHAS = PYALPS(SCALE)
c	     write(6,*) ' ALPAHS (2loop) pythia :', ALPHAS
	   ENDIF
      ENDIF
c      write(6,*) alphas,nf
      RETURN
      END
      FUNCTION ULALPS(Q2)
	Real ULALPS
	Double Precision PYALPS,SCALE
	SCALE=DBLE((Q2))
	ULALPS=REAL(PYALPS(SCALE))
	RETURN
	END
