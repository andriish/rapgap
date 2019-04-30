
      SUBROUTINE LTIMEX(TIME)
C...Interface routine to transfer a call to some machine-dependent
C...routine to get the execution time used since job started.
C...Nice, but not necessary information. Can also be called by user.

      TIME=0.
C...Use of CERN library routine Z007, replace/delete if not available.
ck      CALL TIMEX(TIME)
C-MH activated...
      CALL TIMEX(TIME)
      RETURN
      END
