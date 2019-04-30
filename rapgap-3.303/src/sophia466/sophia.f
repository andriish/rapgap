
C***********************************************************************

       subroutine SOPHIA(L0,P_gam,P_nuc,Imode)
C***********************************************************************
C
C  SOPHIA 1.5
C
C  A. M\"ucke, R. Engel, J.P. Rachen, R.J. Protheroe, and T. Stanev:
C  {\it SOPHIA -- Monte Carlo simulations for photohadronic processes in
C  astrophysics.} BA-99-33 (astro-ph/9903478),
C  submitted to Comp. Phys. Commun.
C
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      SAVE

      DIMENSION P_gam(4),P_nuc(4)

       COMMON /S_RUN/ SQS, S, Q2MIN, XMIN, ZMIN, kb, kt, a1, a2, Nproc
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
       COMMON /S_MASS1/ AM(49), AM2(49)
       COMMON /S_CHP/ S_LIFE(49), ICHP(49), ISTR(49), IBAR(49)
       COMMON /S_CSYDEC/ CBR(102), IDB(49), KDEC(612), LBARP(49)

      CHARACTER NAMPRES*6
      COMMON /RES_PROP/ AMRES(9), SIG0(9),WIDTH(9),
     +                    NAMPRES(0:9)

      CHARACTER NAMPRESp*6
      COMMON /RES_PROPp/ AMRESp(9), BGAMMAp(9),WIDTHp(9),
     +                    RATIOJp(9),NAMPRESp(0:9)

      CHARACTER NAMPRESn*6
      COMMON /RES_PROPn/ AMRESn(9), BGAMMAn(9),WIDTHn(9),
     +                    RATIOJn(9),NAMPRESn(0:9)

       DOUBLE PRECISION P_sum(4),PC(4),GamBet(4)

c       DATA pi /3.141593D0/
       DATA pm /0.93827D0/
       DATA IRESMAX /9/
       DATA Icount / 0 /

c****** INPUT **************************************************
c E0 = energy of incident proton (in lab frame) [in GeV]
c eps = energy of incident photon [in GeV] (in lab frame)
c theta = angle between incident proton and photon [in degrees]
c L0 = code number of the incident nucleon
c***************************************************************
c** calculate eps_prime = photon energy in nuclear rest frame,
c**             sqrt(s) = CMF energy of the N\gamma-system

       Esum  = P_nuc(4)+P_gam(4)
       PXsum = P_nuc(1)+P_gam(1)
       PYsum = P_nuc(2)+P_gam(2)
       PZsum = P_nuc(3)+P_gam(3)
       IQchr = ICHP(1)+ICHP(L0)
       IQbar = IBAR(1)+IBAR(L0)


C  calculate Lorentz boost and rotation
       P_sum(1) = P_nuc(1)+P_gam(1)
       P_sum(2) = P_nuc(2)+P_gam(2)
       P_sum(3) = P_nuc(3)+P_gam(3)
       P_sum(4) = P_nuc(4)+P_gam(4)

       s = (P_sum(4)-P_sum(3))*(P_sum(4)+P_sum(3))
     &     -P_sum(1)**2-P_sum(2)**2

       sqsm = sqrt(s)

       eps_prime = (s-pm*pm)/2.D0/pm
C  Lorentz transformation into c.m. system
      DO I=1,4
        GamBet(I) = P_sum(I)/sqsm
      ENDDO   
C  calculate rotation angles
      IF(GamBet(4).lt.1.d5) then
C  transform nucleon vector
        GamBet(1) = -GamBet(1)
        GamBet(2) = -GamBet(2)
        GamBet(3) = -GamBet(3)
        CALL PO_ALTRA(GamBet(4),GamBet(1),GamBet(2),GamBet(3),
     &                P_nuc(1),P_nuc(2),P_nuc(3),P_nuc(4),Ptot,
     &                PC(1),PC(2),PC(3),PC(4))
        GamBet(1) = -GamBet(1)
        GamBet(2) = -GamBet(2)
        GamBet(3) = -GamBet(3)
C  rotation angle: nucleon moves along +z
        COD = PC(3)/Ptot
        SID = SQRT(PC(1)**2+PC(2)**2)/Ptot
        COF = 1.D0
        SIF = 0.D0
        IF(Ptot*SID.GT.1.D-5) THEN
          COF=PC(1)/(SID*Ptot)
          SIF=PC(2)/(SID*Ptot)
          Anorf=SQRT(COF*COF+SIF*SIF)
          COF=COF/Anorf
          SIF=SIF/Anorf
        ENDIF
      else
        COD = 1.D0
        SID = 0.D0
        COF = 1.D0
        SIF = 0.D0
      endif

c... check for threshold:
       sth = 1.1646D0       
       if (s.lt.sth) then
        print*,'input energy below threshold for photopion production !'
        print*,'sqrt(s) = ',sqrt(s)
        NP = 0
        RETURN
       endif

 200  continue
      Icount = Icount+1
      Imode = 0

c*******************************************************************
c decide which process occurs:                                   ***
c (1) decay of resonance                                         ***
c (2) direct pion production (interaction of photon with         *** 
c     virtual pions in nucleon cloud) and diffractive scattering ***
c (3) multipion production                                       ***
c*******************************************************************

       call dec_inter3(eps_prime,Imode,L0)

c*********************************************
c******* PARTICLE PRODUCTION *****************
c*********************************************
c  42   continue
       if (Imode.le.5) then
c... direct/multipion/diffractive scattering production channel:
        call GAMMA_H(sqsm,L0,Imode,Ifbad)
        if(Ifbad.ne.0) then
          print *,' eventgen: simulation of particle production failed'
          goto 200
        endif
       else if (Imode.eq.6) then
c... Resonances:
c... decide which resonance decays with ID=IRES in list:  
c... IRESMAX = number of considered resonances = 9 so far 
       IRES = 0
 46    call dec_res2(eps_prime,IRES,IRESMAX,L0)
       Nproc = 10+IRES
       call dec_proc2(eps_prime,IPROC,IRANGE,IRES,L0)
c 2-particle decay of resonance in CM system:
       NP = 2
       call res_decay3(IRES,IPROC,IRANGE,s,L0,nbad)
       if (nbad.eq.1) then
         print *,' eventgen: event rejected by res_decay3'
         goto 46
       endif
       call DECSIB
       else
        print*,'invalid Imode !!'
        STOP
       endif

c... consider only stable particles:
 18     istable=0
        do 16 i=1,NP
         if (abs(LLIST(i)).lt.10000) then
          istable = istable+1
          LLIST(istable) = LLIST(i)
          P(istable,1) = P(i,1)
          P(istable,2) = P(i,2)
          P(istable,3) = P(i,3)
          P(istable,4) = P(i,4)
          P(istable,5) = P(i,5)
         endif
  16    continue
        if (NP.gt.istable) then
         do i=istable+1,NP
          LLIST(i) = 0
          P(i,1) = 0.
          P(i,2) = 0.
          P(i,3) = 0.
          P(i,4) = 0.
          P(i,5) = 0.
         enddo
        endif
        NP = istable       

c***********************************************
c transformation from CM-system to lab-system: *
c***********************************************

      DO I=1,NP
        CALL PO_TRANS(P(I,1),P(I,2),P(I,3),COD,SID,COF,SIF,
     &    PC(1),PC(2),PC(3))
        PC(4) = P(I,4)
        CALL PO_ALTRA(GamBet(4),GamBet(1),GamBet(2),GamBet(3),
     &    PC(1),PC(2),PC(3),PC(4),Ptot,
     &    P(I,1),P(I,2),P(I,3),P(I,4))
      ENDDO

      call check_event(Icount,Esum,PXsum,PYsum,PZsum,IQchr,IQbar,Irej)
      if(Irej.ne.0) then
        print *,' eventgen: event rejected by check_event'
        goto 200
      endif

      return

      END
