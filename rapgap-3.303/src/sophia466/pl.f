

      DOUBLE PRECISION function Pl(x,xth,xmax,alpha)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       SAVE

       if (xth.gt.x) then
        Pl = 0.
        RETURN
       endif

       a = alpha*xmax/xth
       prod1 = ((x-xth)/(xmax-xth))**(a-alpha)
       prod2 = (x/xmax)**(-a)
       Pl = prod1*prod2

       END
