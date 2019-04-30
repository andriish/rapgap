************************************************************************
      REAL FUNCTION RKNO (RMASS)

      IMPLICIT NONE
************************************************************************
*
*     Purpose: Determine mean number of pions given by KNO scaling
*     --------
*
*     Input:    RMASS:  Mass available for pion production
*     ------
*
*     Output:   Function value: mean pion number
*     -------
*
*     Author:   B. List  27.1.94
*     -------
*
*     Reference: GOULIANOS, K. (1983): Phys. Rep. 101, 169-219.
*     ---------- Fortmulae 46 to 48
*
************************************************************************
      REAL RMASS,RLNM2 
* Rwarning RLNM2 undefined!

      IF (RMASS .LT. 0.0) THEN
        RKNO = 1.0
      ELSE IF (RMASS .LT. 1.0) THEN
        RKNO = 2.0*SQRT (RMASS)
      ELSE
        RLNM2 = 2.0*LOG (RMASS)
        RKNO = 2.0 + RLNM2*(0.13 + 0.18*RLNM2)
      END IF

      RETURN
      END
