

c*****************************
c*** List of SUBROUTINES *****
C*****************************



      SUBROUTINE INITIAL(L0)

c*******************************************************************
c initialization routine
c*******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      SAVE
      COMMON /RES_PROP/ AMRES(9),SIG0(9),WIDTH(9), 
     +                    NAMPRES(0:9)
      COMMON /RES_PROPp/ AMRESp(9), BGAMMAp(9),WIDTHp(9),  
     +                    RATIOJp(9),NAMPRESp(0:9)
      COMMON /RES_PROPn/ AMRESn(9), BGAMMAn(9),WIDTHn(9),  
     +                    RATIOJn(9),NAMPRESn(0:9)
      COMMON /S_MASS1/ AM(49), AM2(49)
      CHARACTER NAMPRESp*6, NAMPRESn*6
      CHARACTER NAMPRES*6
C...for interface to DJANGOH
      COMMON /SPPASS/ NSOPH,NSPOUT,NFAILP,NSPACC

      NSOPH=0
      NSPOUT=0
      NFAILP=0

*      print *,'settings of resonance parameters:'
*      print *,
*    &   'name, nom.mass [GeV], width [GeV], sigma0 [mub], bgamma*10^3'

       if (L0.eq.13) then
       do i=1,9
        SIG0(i) = 4.893089117D0/AM2(13)*RATIOJp(i)*BGAMMAp(i)
        AMRES(i) = AMRESp(i)
        WIDTH(i) = WIDTHp(i)
        NAMPRES(i) = NAMPRESp(i)
*      print*,NAMPRES(i), real(AMRES(i)), real(WIDTH(i)),
*    &          real(SIG0(i)), real(BGAMMAp(i))
       enddo
       endif

       if (L0.eq.14) then
       do i=1,9
        SIG0(i) = 4.893089117D0/AM2(14)*RATIOJn(i)*BGAMMAn(i)
        AMRES(i) = AMRESn(i)
        WIDTH(i) = WIDTHn(i)
        NAMPRES(i) = NAMPRESn(i)
*      print*,NAMPRES(i), real(AMRES(i)), real(WIDTH(i)),
*    &         real(SIG0(i)), real(BGAMMAn(i))
       enddo
       endif

       RETURN
       END
