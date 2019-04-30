
      DOUBLE PRECISION function probangle(IRES,L0,z)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       SAVE

c********************************************************************
c probability distribution for scattering angle of given resonance **
c IRES and incident nucleon L0 ;                                   **
c z is cosine of scattering angle in CMF frame                     **
c********************************************************************

       if (IRES.eq.4.or.IRES.eq.5.or.IRES.eq.2) then  
c ... N1535 andf N1650 decay isotropically. 
        probangle = 0.5D0 
        return
       endif

       if (IRES.eq.1) then
c ... for D1232:  
        probangle =  0.636263D0 - 0.408790D0*z*z
        return
       endif

       if (IRES.eq.3.and.L0.eq.14) then
c ... for N1520 and incident n: 
        probangle =  0.673669D0 - 0.521007D0*z*z
        return
       endif

       if (IRES.eq.3.and.L0.eq.13) then
c ... for N1520 and incident p: 
        probangle =  0.739763D0 - 0.719288D0*z*z
        return
       endif

       if (IRES.eq.6.and.L0.eq.14) then
c ... for N1680 (more precisely: N1675) and incident n: 
        q=z*z
        probangle = 0.254005D0 + 1.427918D0*q - 1.149888D0*q*q
        return
       endif


       if (IRES.eq.6.and.L0.eq.13) then
c ... for N1680 and incident p: 
        q=z*z
        probangle = 0.189855D0 + 2.582610D0*q - 2.753625D0*q*q
        return
       endif

      if (IRES.eq.7) then
c ... for D1700:  
       probangle =  0.450238D0 + 0.149285D0*z*z
       return
      endif


      if (IRES.eq.8) then
c ... for D1905:  
       q=z*z
       probangle = 0.230034D0 + 1.859396D0*q - 1.749161D0*q*q
       return
      endif


      if (IRES.eq.9) then
c ... for D1950:  
       q=z*z
       probangle = 0.397430D0 - 1.498240D0*q + 5.880814D0*q*q
     &                - 4.019252D0*q*q*q
       return
      endif

      print*,'error in function probangle !'
      STOP
      END
