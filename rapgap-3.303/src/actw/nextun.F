C================================================================

      FUNCTION NEXTUN()
C                                    Returns an unallocated FORTRAN i/o unit.
      LOGICAL EX
C
      DO 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         IF (.NOT. EX) then
            nextun = n
            RETURN
         end if
   10 CONTINUE
      RETURN
C               *************************
      END
