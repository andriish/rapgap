      FUNCTION RLU()
      Double Precision draprn
      EXTERNAL draprn
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      IF(FIRST) THEN
         write(6,*) 'call rlu = draprn'
         FIRST = .FALSE.
      ENDIF
      RLU =  draprn()
      RETURN
      END

