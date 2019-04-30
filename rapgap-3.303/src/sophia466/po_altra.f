

      SUBROUTINE PO_ALTRA(GA,BGX,BGY,BGZ,PCX,PCY,PCZ,EC,P,PX,PY,PZ,E)
C*********************************************************************
C
C    arbitrary Lorentz transformation
C
C    (taken from PHOJET v1.12, R.E. 08/98)
C
C*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      EP=PCX*BGX+PCY*BGY+PCZ*BGZ
      PE=EP/(GA+1.D0)+EC
      PX=PCX+BGX*PE
      PY=PCY+BGY*PE
      PZ=PCZ+BGZ*PE
      P=SQRT(PX*PX+PY*PY+PZ*PZ)
      E=GA*EC+EP

      END
