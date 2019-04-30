

      DOUBLE PRECISION function Ef(x,th,w)

      IMPLICIT DOUBLE PRECISION (A-M,O-Z)
      IMPLICIT INTEGER (N)

       SAVE

       wth = w+th
       if (x.le.th) then
        Ef = 0.
        RETURN
       else if (x.gt.th.and.x.lt.wth) then
        Ef = (x-th)/w
        RETURN
       else if (x.ge.wth) then
        Ef = 1.
        RETURN
       else
        print*,'error in function EF'
        STOP
       endif

       END
