      SUBROUTINE HS46INI
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C...EXTERNAL statement links HSBLKD on most machines.
      EXTERNAL hsblkd

      write(6,*) ' HERACLES 4.6.6 BLOCK DATA initialized '
      RETURN
      END
#include "hsblkd.fsave"
