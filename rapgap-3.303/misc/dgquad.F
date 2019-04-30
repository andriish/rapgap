*
* $Id: dgquad.F,v 1.1.1.1 1996/04/01 15:02:15 mclareni Exp $
*
* $Log: dgquad.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:15  mclareni
* Mathlib gen
*
*
*         Double Prec. version only on 32 bit machines.
      FUNCTION DGQUAD(F,A,B,N)
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION X(1),W(1)
      EXTERNAL F

      CALL D107D1(1,F,A,B,N,X,W)
      DGQUAD=X(1)
      RETURN
      END
