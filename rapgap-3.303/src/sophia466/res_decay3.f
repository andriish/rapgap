

       SUBROUTINE RES_DECAY3(IRES,IPROC,IRANGE,s,L0,nbad)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

       SAVE

       COMMON /S_RESp/ CBRRES1p(18),CBRRES2p(36),CBRRES3p(26),
     +  RESLIMp(36),ELIMITSp(9),KDECRES1p(90),KDECRES2p(180),
     +  KDECRES3p(130),IDBRES1p(9),IDBRES2p(9),IDBRES3p(9) 
       COMMON /S_RESn/ CBRRES1n(18),CBRRES2n(36),CBRRES3n(22),
     +  RESLIMn(36),ELIMITSn(9),KDECRES1n(90),KDECRES2n(180),
     +  KDECRES3n(110),IDBRES1n(9),IDBRES2n(9),IDBRES3n(9) 
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
c       COMMON /S_CNAM/ NAMP (0:49)
c      CHARACTER NAMP*6, NAMPRESp*6, NAMPRESn*6

*      external scatangle, proc_twopart

c********************************************************
c  RESONANCE AMD with code number IRES  INTO  M1 + M2
C  PROTON ENERGY E0 [in GeV] IN DMM [in GeV]
C  E1,E2 [in GeV] are energies of decay products
c  LA,LB are code numbers of decay products
c  P(1,1:5),P(2,1:5) are 5-momenta of particles LA,LB;
c  resulting momenta are calculated in CM frame;
c  ANGLESCAT is cos of scattering angle in CM frame
c********************************************************
c** Date: 20/01/98   **
c** correct.:28/04/98**
c** author: A.Muecke **
c**********************

c... determine decay products LA, LB:
        NP = 2
        if (L0.eq.13) then
c ... proton is incident nucleon:
        if (IRANGE.eq.1) then
         LA = KDECRES1p(5*(IPROC-1)+3)
         LB = KDECRES1p(5*(IPROC-1)+4)
        else if (IRANGE.eq.2) then
         LA = KDECRES2p(5*(IPROC-1)+3)
         LB = KDECRES2p(5*(IPROC-1)+4)
        else if (IRANGE.eq.3) then
         LA = KDECRES3p(5*(IPROC-1)+3)
         LB = KDECRES3p(5*(IPROC-1)+4)
        else 
          print*,'error in res_decay3'
        endif
        else if (L0.eq.14) then
c ... neutron is incident nucleon:
        if (IRANGE.eq.1) then
         LA = KDECRES1n(5*(IPROC-1)+3)
         LB = KDECRES1n(5*(IPROC-1)+4)
        else if (IRANGE.eq.2) then
         LA = KDECRES2n(5*(IPROC-1)+3)
         LB = KDECRES2n(5*(IPROC-1)+4)
        else if (IRANGE.eq.3) then
         LA = KDECRES3n(5*(IPROC-1)+3)
         LB = KDECRES3n(5*(IPROC-1)+4)
        else 
          print*,'error in res_decay3'
        endif

        else
         print*,'no valid L0 in RES_DECAY'
         STOP
        endif

        LLIST(1) = LA
        LLIST(2) = LB

c... sample scattering angle:
       call scatangle(anglescat,IRES,L0)
       
c ... 2-particle decay:
        call proc_twopart(LA,LB,sqrt(s),LLIST,P,anglescat,nbad)

        RETURN

        END
