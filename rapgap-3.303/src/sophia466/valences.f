

      SUBROUTINE valences(ip,ival1,ival2)
C**********************************************************************
C
C    valence quark composition of various particles  (R.E. 03/98)
C    (with special treatment of photons)
C
C**********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      SAVE

      if(ip.eq.1) then
c... photon couples to quark:
        if(rndm(0).gt.0.2D0) then
c... photon couples to meson: pion-cloud model
*       if(rndm(0).gt.0.7D0) then
          ival1 = 1
          ival2 = -1
        else
          ival1 = 2
          ival2 = -2
        endif
      else if(ip.eq.6) then
        if(rndm(0).gt.0.5D0) then
          ival1 = 1
          ival2 = -1
        else
          ival1 = 2
          ival2 = -2
        endif
      else if(ip.eq.7) then
        ival1 = 1
        ival2 = -2
      else if(ip.eq.8) then
        ival1 = 2
        ival2 = -1
      else if(ip.eq.13) then
        Xi = rndm(0)
        if(Xi.lt.0.3333D0) then
          ival1 = 12
          ival2 = 1
        else if(Xi.lt.0.6666D0) then
          ival1 = 21
          ival2 = 1
        else
          ival1 = 11
          ival2 = 2
        endif
      else if(ip.eq.14) then
        Xi = rndm(0)
        if(Xi.lt.0.3333D0) then
          ival1 = 12
          ival2 = 2
        else if(Xi.lt.0.6666D0) then
          ival1 = 21
          ival2 = 2
        else
          ival1 = 22
          ival2 = 1
        endif
      endif

      if((ip.lt.13).and.(rndm(0).lt.0.5D0)) then
        k = ival1
        ival1 = ival2
        ival2 = k
      endif

      END
