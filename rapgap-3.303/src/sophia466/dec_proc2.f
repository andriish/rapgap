

      subroutine dec_proc2(x,IPROC,IRANGE,IRES,L0)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       SAVE

c**********************************************************************
c*** decide which decay with ID=IPROC of resonance IRES takes place ***
c**********************************************************************
c** Date: 20/01/98   **
c** correct.: 27/04/98*
c** author: A.Muecke **
c**********************

       COMMON /S_RESp/ CBRRES1p(18),CBRRES2p(36),CBRRES3p(26),
     +  RESLIMp(36),ELIMITSp(9),KDECRES1p(90),KDECRES2p(180),
     +  KDECRES3p(130),IDBRES1p(9),IDBRES2p(9),IDBRES3p(9)
       COMMON /S_RESn/ CBRRES1n(18),CBRRES2n(36),CBRRES3n(22),
     +  RESLIMn(36),ELIMITSn(9),KDECRES1n(90),KDECRES2n(180),
     +  KDECRES3n(110),IDBRES1n(9),IDBRES2n(9),IDBRES3n(9)
       DIMENSION prob_sum(0:9)

c      x = eps_prime
c ... choose arrays /S_RESp/ for charged resonances,
c ...        arrays /S_RESn/ for neutral resonances
       if (L0.eq.13) then
c ... charged resonances:

       r = RNDM(0)
c... determine the energy range of the resonance:
       nlim = ELIMITSp(IRES)
       istart = (IRES-1)*4+1
       if (nlim.gt.0) then
         do ie=istart,nlim-2+istart
           reslimp1 = RESLIMp(ie)
           reslimp2 = RESLIMp(ie+1)
          if (x.le.reslimp2.and.x.gt.reslimp1) then
           IRANGE = ie+1-istart
          endif
         enddo
       else
         irange = 1
  13   endif



       IPROC = -1
       i = 0
       prob_sum(0) = 0.D0

       if (IRANGE.eq.1) then
        j = IDBRES1p(IRES)-1
        if (j.eq.-1) then
         print*,'invalid resonance in energy range 1'
        endif
 10     continue
        j = j+1
        i = i+1
        prob_sum(i) = CBRRES1p(j)
        if (r.ge.prob_sum(i-1).and.r.lt.prob_sum(i)) then
         IPROC = j
        endif
        if (prob_sum(i).lt.1.D0) goto 10
        if (r.eq.1.D0) IPROC = j
        if (IPROC.eq.-1) then
         print*,'no resonance decay possible !'
        endif

       else if (IRANGE.eq.2) then
        j = IDBRES2p(IRES)-1
        if (j.eq.-1) then
         print*,'invalid resonance in energy range 2'
        endif
 11     continue
        j = j+1
        i = i+1
        prob_sum(i) = CBRRES2p(j)
        if (r.ge.prob_sum(i-1).and.r.lt.prob_sum(i)) then
         IPROC = j
        endif
        if (prob_sum(i).lt.1.D0) goto 11
        if (r.eq.1.D0) IPROC = j
        if (IPROC.eq.-1) then
         print*,'no resonance decay possible !'
        endif

       else if (IRANGE.eq.3) then
        j = IDBRES3p(IRES)-1
        if (j.eq.-1) then
         print*,'invalid resonance in energy range 3'
        endif
 12     continue
        j = j+1
        i = i+1
        prob_sum(i) = CBRRES3p(j)
        if (r.ge.prob_sum(i-1).and.r.lt.prob_sum(i)) then
         IPROC = j
        endif
        if (prob_sum(i).lt.1.D0) goto 12
        if (r.eq.1.D0) IPROC = j
        if (IPROC.eq.-1) then
         print*,'no resonance decay possible !'
        endif

        else
         print*,'invalid IRANGE in DEC_PROC2'
        endif

       RETURN


         else if (L0.eq.14) then
c ... neutral resonances:

       r = RNDM(0)
c... determine the energy range of the resonance:
       nlim = ELIMITSn(IRES)
       istart = (IRES-1)*4+1
       if (nlim.gt.0) then
         do ie=istart,nlim-2+istart
          if (x.le.RESLIMn(ie+1).and.x.gt.RESLIMn(ie)) then
           IRANGE = ie+1-istart
          endif
         enddo
       else
         irange = 1
       endif


       IPROC = -1
       i = 0
       prob_sum(0) = 0.D0

       if (IRANGE.eq.1) then
        j = IDBRES1n(IRES)-1
        if (j.eq.-1) then
         print*,'invalid resonance in this energy range'
        endif
 20     continue
        j = j+1
        i = i+1
        prob_sum(i) = CBRRES1n(j)
        if (r.ge.prob_sum(i-1).and.r.lt.prob_sum(i)) then
         IPROC = j
        endif
        if (prob_sum(i).lt.1.D0) goto 20
        if (r.eq.1.D0) IPROC = j
        if (IPROC.eq.-1) then
         print*,'no resonance decay possible !'
        endif

       else if (IRANGE.eq.2) then
        j = IDBRES2n(IRES)-1
        if (j.eq.-1) then
         print*,'invalid resonance in this energy range'
        endif
 21     continue
        j = j+1
        i = i+1
        prob_sum(i) = CBRRES2n(j)
        if (r.ge.prob_sum(i-1).and.r.lt.prob_sum(i)) then
         IPROC = j
        endif
        if (prob_sum(i).lt.1.D0) goto 21
        if (r.eq.1.) IPROC = j
        if (IPROC.eq.-1) then
         print*,'no resonance decay possible !'
        endif

       else if (IRANGE.eq.3) then
        j = IDBRES3n(IRES)-1
        if (j.eq.-1) then
         print*,'invalid resonance in this energy range'
        endif
 22     continue
        j = j+1
        i = i+1
        prob_sum(i) = CBRRES3n(j)
        if (r.ge.prob_sum(i-1).and.r.lt.prob_sum(i)) then
         IPROC = j
        endif
        if (prob_sum(i).lt.1.D0) goto 22
        if (r.eq.1.D0) IPROC = j
        if (IPROC.eq.-1) then
         print*,'no resonance decay possible !'
        endif

        else
         print*,'invalid IRANGE in DEC_PROC2'
        endif

       RETURN

       else
        print*,'no valid L0 in DEC_PROC !'
        STOP
       endif

       END
