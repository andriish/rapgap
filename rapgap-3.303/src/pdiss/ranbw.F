************************************************************************
      REAL FUNCTION RANBW (ER, GAMMA, EMIN, EMAX)
      
      IMPLICIT NONE
************************************************************************
*
*     Purpose:   Generate random number with BREIT-WIGNER distribution
*     --------
*
*     Input:     ER:     Maximum of distribution
*     ------     GAMMA:  Width of distribution
*                EMIN:   Minimal value of RANBW
*                EMAX:   Maximal value of RANBW
*
*     Output:    Function value: Random number between EMIN and EMAX
*     -------    with BREIT-WIGNER distribution:
*                1 / ((E - ER)**2 + GAMMA**2/4)
*
*     Author:    B. List,  18.7.93
*     -------
*
************************************************************************
      REAL*8 ER, GAMMA, EMIN, EMAX,A,B,E,DRAPRN
      

      IF (GAMMA .LT. 1.0E-3*ER) THEN
         RANBW = SNGL(ER)
      ELSE
         A = ATAN (2.0*(EMAX - ER)/GAMMA)
         B = ATAN (2.0*(EMIN - ER)/GAMMA)
         E = ER + 0.5*GAMMA*TAN (DRAPRN ()*(A - B) + B)
         IF (E .LE. EMAX) THEN
            RANBW = SNGL(E)
         ELSE
            RANBW = SNGL(EMAX)
         END IF
      END IF

      RETURN
      END
