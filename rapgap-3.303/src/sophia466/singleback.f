
c***********************************************************
C calculates functions for crossection of direct channel 
c NOT isospin-corrected, simply a samelsurium of functions
c x is eps_prime in GeV (see main program)
C (see thesis of J.Rachen, p.45ff)
c note: neglect strange- and eta-channel
C***********************************************************
c** Date: 27/04/98   **
c** last chg:23/05/98**
c** author: A.Muecke **
c**********************
c

       DOUBLE PRECISION FUNCTION singleback(x)
c****************************
c SINGLE PION CHANNEL
c****************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       SAVE
       singleback = 92.7D0*Pl(x,.152D0,.25D0,2.D0)

       END
