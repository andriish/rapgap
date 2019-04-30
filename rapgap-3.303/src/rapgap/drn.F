c get draprn for bases/spring
      FUNCTION DRN(ISEED)
	Implicit None
	Integer ISEED
      double precision drn1,DRN
      Double Precision draprn
	EXTERNAL draprn
      DRN1= draprn()
      DRN = DRN1
	RETURN
	END
