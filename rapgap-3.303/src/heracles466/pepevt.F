      SUBROUTINE PEPEVT
C**********************
C
C     print out the EPEVT common in a suitable form
C
C     Peter Lanius
C     07/02/92
C
C*********************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
#include "hepevt.inc"
C---EVENT CHARACTERISTICS IN COMMON EPEVT
C   ELECTRON:      IDHEP(1), PHEP(I,1)
C   QUARK:         IDHEP(2), PHEP(I,2)
C   PHOTON:        IDHEP(3), PHEP(I,3)

      IF( NHEP .EQ. 0 ) THEN
         WRITE(*,*) '+++++++++++ PEPEVT: NHEP = 0, NO OUTPUT ++++++++++'
     +   //'+'
      END IF

      IF( NHEP .GT. 0 ) THEN
         WRITE(6,10000)
         DO 10  I = 1, NHEP
            WRITE(6,10100) I, IDHEP(I), ISTHEP(I), JMOHEP(1,I), PHEP(1,
     +      I), PHEP(2,I), PHEP(3,I), PHEP(4,I)
   10    CONTINUE
      END IF

10000 FORMAT(//,120('-'),/,6X,'Number',2X,'PDG-code',3X,
     +       'St-code',1X,'1.Mother',7X,'PX',13X,'PY',
     +       13X,'PZ',10X,'Energy',/,/,120('-'),/)
10100 FORMAT(4I10,4D15.4)

      RETURN
      END
