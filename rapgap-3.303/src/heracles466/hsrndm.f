      DOUBLE PRECISION FUNCTION HSRNDM()
      Double Precision draprn
      EXTERNAL draprn
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      IF(FIRST) THEN
         write(6,*) ' change of random number generator:               '
     +   //'              call hsrndm = draprn'
         FIRST = .FALSE.
      ENDIF
      HSRNDM = draprn()
      RETURN
      END
