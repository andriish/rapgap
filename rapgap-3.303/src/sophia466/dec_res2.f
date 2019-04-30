

      subroutine dec_res2(eps_prime,IRES,IRESMAX,L0)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       SAVE

c*****************************************************************************
c*** decides which resonance with ID=IRES in list takes place at eps_prime ***
c*****************************************************************************
c** Date: 20/01/98   **
c** author: A.Muecke **
c**********************

       DIMENSION prob_sum(9)
c       CHARACTER NAMPRESp*6,NAMPRESn*6


c*** sum of all resonances:
       sumres = 0.D0
       do 12 j=1,IRESMAX
        j10 = j+10
        sumres = sumres+crossection(eps_prime,j10,L0)
        prob_sum(j) = sumres
  12   continue


       r = RNDM(0)

       IRES = 0
       i = 0
       prob = 0.D0
 10    continue
       i = i+1
       probold = prob
       prob = prob_sum(i)/sumres
       if (r.ge.probold.and.r.lt.prob) then
         IRES = i
         RETURN
       endif
       if (i.lt.IRESMAX) goto 10
       if (r.eq.1.D0) IRES = i
       if (IRES.eq.0) then
         print*,'no resonance possible !'
         STOP
       endif

       RETURN

       END
