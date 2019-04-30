      FUNCTION RNDM()
	Implicit None
	Double Precision RNDM
      Double Precision draprn
      EXTERNAL draprn
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      IF(FIRST) THEN
         write(6,*) 'call rndm = draprn'
         FIRST = .FALSE.
      ENDIF
      RNDM = draprn()
      RETURN
      END
