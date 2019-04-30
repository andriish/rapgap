

      subroutine scatangle(anglescat,IRES,L0)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       SAVE

c*******************************************************************
c This routine samples the cos of the scattering angle for a given *
c resonance IRES and incident nucleon L0; it is exact for         **
c one-pion decay channel and if there is no                       **
c other contribution to the cross section from another resonance  **
c and an approximation for an overlay of resonances;              **
c for decay channels other than the one-pion decay a isotropic    **
c distribution is used                                            **
c*******************************************************************
c** Date: 16/02/98   **
c** author: A.Muecke **
c**********************

       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb

c ... use rejection method for sampling:
       LA = LLIST(1)
       LB = LLIST(2)
  10   continue
       r = RNDM(0)
c*** sample anglescat random between -1 ... 1 **
      anglescat = 2.D0*(r-0.5D0) 
c ... distribution is isotropic for other than one-pion decays:
       if ((LA.eq.13.or.LA.eq.14).and.LB.ge.6.and.LB.le.8) then
        prob = probangle(IRES,L0,anglescat)
       else
        prob = 0.5D0
       endif
       r = RNDM(0)
       if (r.le.prob) then
          RETURN
        else
         goto 10
       endif       
 12   continue

       END
