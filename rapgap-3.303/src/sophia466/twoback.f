

       DOUBLE PRECISION FUNCTION twoback(x)
c*****************************
c TWO PION PRODUCTION
c*****************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       SAVE
       twoback = 37.7D0*Pl(x,.4D0,.6D0,2.D0)

       END
