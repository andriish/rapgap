

      SUBROUTINE PROC_TWOPART(LA,LB,AMD,Lres,Pres,costheta,nbad)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      COMMON /S_MASS1/ AM(49), AM2(49)
      COMMON /RES_FLAG/ FRES(49),XLIMRES(49)
      SAVE
      DIMENSION Pres(2000,5),Lres(2000)

c***********************************************************
c  2-particle decay of CMF mass AMD INTO  M1 + M2
C  NUCLEON ENERGY E0 [in GeV];
C  E1,E2 [in GeV] are energies of decay products
c  LA,LB are code numbers of decay products
c  P1(1:5),P2(1:5) are 5-momenta of particles LA,LB;
c  resulting momenta are calculated in CM frame;
c  costheta is cos of scattering angle in CM frame
c  this program also checks if the resulting particles are
c  resonances; if yes, it is also allowed to decay a
c  mass AMD < M1 + M2 by using the width of the resonance(s)
c***********************************************************
c** Date: 20/01/98   **
c** correct.:19/02/98**
c** author: A.Muecke **
c**********************

        nbad = 0
c        ND = 2
        SM1 = AM(LA)
        if (LB.eq.0) then
c... for Landau-Model:
         SM2 = 2.D0*AM(7)
        else
         SM2 = AM(LB)
        endif
	E1 = (AMD*AMD + SM1*SM1 - SM2*SM2)/AMD/2.D0
	E2 = (AMD*AMD + SM2*SM2 - SM1*SM1)/AMD/2.D0
c... check if SM1+SM2 < AMD:
        if ((SM1+SM2).gt.AMD) then
c... if one of the decay products is a resonance, this 'problem' can
c    be solved by using a reduced mass for the resonance and assume that
c    this resonance is produced at its threshold;
         if (FRES(LA).eq.1.D0) then
c ...      particle LA is a resonance:
          SM1 = AMD-SM2
	  E1 = SM1
	  E2 = AMD-E1
         if (E1.lt.XLIMRES(LA).or.E2.lt.XLIMRES(LB)) nbad = 1
         endif
        if (FRES(LB).eq.1.D0) then
c ...      particle LB is a resonance:
          SM2 = AMD-SM1
	  E2 = SM2
         E1 = AMD-E2
          if (E1.lt.XLIMRES(LA).or.E2.lt.XLIMRES(LB)) nbad = 1
         endif
c ...     both particles are NOT resonances: -> error !  
         if (FRES(LA).eq.0.D0.and.FRES(LB).eq.0.D0) then
          print*,'SM1 + SM2 > AMD in PROC_TWOPART',SM1,SM2,AMD,LA,LB
          STOP
         endif
        endif

       if (nbad.eq.0) then
	PC = SQRT((E1*E1 - SM1*SM1))
        Pres(1,4) = E1
        Pres(2,4) = E2
        Pres(1,5) = SM1
        Pres(2,5) = SM2
        
        
C *********************************************************
c theta is scattering angle in CM frame: 
        r = RNDM(0)
        P1Z= PC*costheta
        P2Z=-PC*costheta

        P1X = sqrt(r*(PC*PC-P1Z*P1Z))
        P2X = sqrt(r*(PC*PC-P2Z*P2Z))
        P1Y = sqrt((1.D0-r)*(PC*PC-P1Z*P1Z))
        P2Y = sqrt((1.D0-r)*(PC*PC-P2Z*P2Z))
        if(RNDM(0).lt.0.5D0) then
          P1X = -P1X
        else
          P2X = -P2X
        endif
        if(RNDM(0).lt.0.5D0) then
          P1Y = -P1Y
        else
          P2Y = -P2Y
        endif

        Pres(1,1) = P1X
        Pres(1,2) = P1Y
        Pres(1,3) = P1Z
        Pres(2,1) = P2X
        Pres(2,2) = P2Y
        Pres(2,3) = P2Z
        Lres(1) = LA
        Lres(2) = LB
       endif

        RETURN
 
        END
