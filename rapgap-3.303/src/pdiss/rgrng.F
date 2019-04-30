*-- Author :    Alexander Fedotov  

************************************************************************
      FUNCTION RGRNG ()
************************************************************************
************************************************************************
*
*
*  Purpose: Generate a number distributed according to a Gaussian
*  -------  with mean = 0 and sigma = 1 (the normal distribution)
*  Author:  A.Fedotov  
************************************************************************

      PARAMETER ( PI = 3.141 592 7 , TWOPI = 2. * PI)
      Double Precision U(2)
      DATA I /0/
*  Initialise the variable to calm down a compiler:
      DATA Z2 /0./
*--------------------

      IF (I .EQ. 0) THEN

*  Generate two gauss. random numbers, Z1 and Z2, and return Z1

         CALL draprnv (U, 2)
         SQ = SQRT ( -2. * DLOG (U(2)))
         U12PI = TWOPI * U(1)
         Z1 = SIN (U12PI) * SQ
         Z2 = COS (U12PI) * SQ

         RGRNG = Z1
         I = 1

      ELSE

*  Return Z2 , generated at previous call

         RGRNG = Z2
         I = 0

      ENDIF

      RETURN

      RETURN
      END

