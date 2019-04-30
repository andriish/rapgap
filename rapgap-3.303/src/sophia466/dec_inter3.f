


      subroutine dec_inter3(eps_prime,Imode,L0)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      SAVE

       DOUBLE PRECISION RNDM
       external RNDM

c*** decides which process takes place at eps_prime ********
c (6) excitation/decay of resonance                      ***
c (2) direct pion production: N\gamma --> N \pi          *** 
c (3) direct pion production: N\gamma --> \Delta \pi     *** 
c (1) diffractive scattering: N\gamma --> N \rho         ***
c (4) diffractive scattering: N\gamma --> N \omega       ***
c (0) multipion production (fragmentation)               ***
c (5) fragmentation in resonance region                  ***
c***********************************************************
c** Date: 15/04/98   **
c** author: A.Muecke **
c**********************
       tot = crossection(eps_prime,3,L0)
       if (tot.eq.0.) tot = 1.D0
       prob1 = crossection(eps_prime,1,L0)/tot
       prob2 = crossection(eps_prime,7,L0)/tot
       prob3 = crossection(eps_prime,2,L0)/tot
       prob4 = crossection(eps_prime,8,L0)/tot
       prob5 = crossection(eps_prime,9,L0)/tot
       prob6 = crossection(eps_prime,0,L0)/tot
       prob7 = 1.D0
       rn = RNDM(0)
       if (rn.lt.prob1) then
        Imode = 6
c ... --> resonance decay
        RETURN
       else if (prob1.le.rn.and.rn.lt.prob2) then
        Imode = 2
c ... --> direct channel: N\gamma --> N\pi
        RETURN
       else if (prob2.le.rn.and.rn.lt.prob3) then
        Imode = 3
c ... --> direct channel: N\gamma --> \Delta \pi
        RETURN
       else if (prob3.le.rn.and.rn.lt.prob4) then
        Imode = 1
c ... --> diffractive scattering: N\gamma --> N \rho
        RETURN
       else if (prob4.le.rn.and.rn.lt.prob5) then
        Imode = 4
c ... --> diffractive scattering: N\gamma --> N \omega
        RETURN
       else if (prob5.le.rn.and.rn.lt.prob6) then
        Imode = 5
c ... --> fragmentation (2) in resonance region
        return
       else if (prob6.le.rn.and.rn.lt.1.) then
        Imode = 0
c ... --> fragmentation mode/multipion production
        RETURN
       else if (rn.eq.1.D0) then
        Imode = 0
        RETURN
       else
        print*,'error in dec_inter.f !'
        STOP
       endif

        END
