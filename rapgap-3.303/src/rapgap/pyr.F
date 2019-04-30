*CMZ :  1.00/07 05/10/2001  10.50.57  by  Hannes Jung
*CMZ :  1.00/04 22/09/2001  20.15.50  by  Hannes Jung
*CMZ :  1.00/01 12/06/2001  15.00.22  by  Hannes Jung
      FUNCTION PYR(IDUMMY)

      double precision pyr,draprn
	external draprn
	
      PYR= draprn(idummy)

      RETURN
      END
