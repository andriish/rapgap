      SUBROUTINE LSCALE(INFIN,QMAX)
      IMPLICIT NONE
      Double Precision QMAX
#include "rgfull.inc"
#include "rgpara1.inc"
#include "rgparas.inc"
#include "rglujets.inc"
#include "rgdisdif.inc"
#include "rgnlo.inc"
      INTEGER INFIN
      IF(IFULL.EQ.0) THEN
         IF(IPRO.EQ.12) THEN
            QMAX = Q2
         ELSE
            QMAX = P(27,5)**2
            IF(IPRO.EQ.18) THEN
               QMAX = Q2Q
                ENDIF
            IF(INFIN.LT.0) THEN
               IF(IPRO.EQ.18) THEN
                  QMAX = Q2Q
               ELSE
                  QMAX=DMAX1(DABS(-Q2-2.D0*DOT1(25,21)), DABS(-Q2-
     +            2.D0* DOT1(26,21)))
               ENDIF
            ENDIF
         ENDIF
      ELSE
         IF(IPRO.EQ.12) THEN
	      IF(IBS.EQ.0) THEN
               QMAX = 4.d0*PT2CUT(13)
		ELSE
		   QMAX = Q2
		ENDIF
         ELSE
            QMAX = P(27,5)**2
            IF(IPRO.EQ.99) QMAX = P(28,5)**2
            IF(INFIN.LT.0) QMAX=DMAX1(DABS(-Q2-2.D0*DOT1(25,21)),
     +      DABS(-Q2-2.D0*DOT1(26,21)))
            IF(IPRO.EQ.99.and.INFIN.LT.0) THEN
               QMAX=DMAX1(DABS(-Q2-2.D0*DOT1(25,21)), DABS(-Q2-2.D0*
     +         DOT1(27,21)))
            ENDIF
         ENDIF
      ENDIF
      QMAX=SQRT(QMAX)
c	write(6,*) ' lscale ',qmax,infin,ifull
      RETURN
      END
