

      DOUBLE PRECISION function crossection(x,NDIR,NL0)

      IMPLICIT DOUBLE PRECISION (A-M,O-Z)
      IMPLICIT INTEGER (N)

      SAVE

      CHARACTER NAMPRES*6
      COMMON /RES_PROP/ AMRES(9), SIG0(9),WIDTH(9), 
     +                    NAMPRES(0:9)

      DIMENSION sig_res(9)

       external breitwigner, Ef, singleback, twoback

       DATA pm /0.93827D0/
       DATA sth /1.1646D0/

c*****************************************************
C calculates crossection of N-gamma-interaction
C (see thesis of J.Rachen, p.45ff and corrections 
C  report from 27/04/98, 5/05/98, 22/05/98 of J.Rachen)
C*****************************************************
c** Date: 20/01/98   **
c** correct.:27/04/98**
c** update: 23/05/98 **
c** author: A.Muecke **
c**********************
c
c x = eps_prime in GeV
c       x = eps_prime
       s = pm*pm+2.D0*pm*x
       
       if (s.lt.sth) then
c       if (x.lt.0.15275) then
        crossection = 0.
c        print*,'BELOW P\GAMMA THRESHOLD !'
        RETURN
       endif
       if (x.gt.10.D0) then
c only multipion production:
        cross_res = 0.D0
        cross_dir = 0.D0
        cross_dir1 = 0.D0
        cross_dir2 = 0.D0
        goto 10
       endif

c****************************
c RESONANCES:
c****************************  

      cross_res = 0.D0

       cross_res = breitwigner(SIG0(1),WIDTH(1),AMRES(1),x)
     &              *Ef(x,0.152D0,0.17D0)
       sig_res(1) = cross_res
      DO N=2,9

        sig_res(N) = breitwigner(SIG0(N),WIDTH(N),AMRES(N),x)
     &              *Ef(x,0.15D0,0.38D0)
        cross_res = cross_res + sig_res(N)

      ENDDO

c****************************
c DIRECT CHANNEL:
c****************************  

       if((x.gt.0.1D0).and.(x.lt.0.6D0)) then
         cross_dir1 = singleback(x)
     &               + 40.D0*exp(-(x-0.29D0)**2/0.002D0)
     &               - 15.D0*exp(-(x-0.37D0)**2/0.002D0)
       else
         cross_dir1 = singleback(x)
       endif
       cross_dir2 = twoback(x)

       cross_dir = cross_dir1 + cross_dir2

c****************************
c FRAGMENTATION 2:
c**************************** 
 10   continue 
       if (NL0.eq.13) then
        cross_frag2 = 80.3D0*Ef(x,0.5D0,0.1D0)*(s**(-0.34D0)) 
       else if (NL0.eq.14) then
c        cross_frag2 = 72.3*Ef(x,0.5D0,0.1D0)*(s**(-0.34)) 
        cross_frag2 = 60.2D0*Ef(x,0.5D0,0.1D0)*(s**(-0.34D0))
       endif

c****************************************************
c MULTIPION PRODUCTION/FRAGMENTATION 1 CROSS SECTION
c****************************************************
       if (x.gt.0.85D0) then
        ss1 = (x-.85D0)/.69D0
        if (NL0.eq.13) then
         ss2 = 29.3D0*(s**(-.34D0))+59.3D0*(s**.095D0)
        else if (NL0.eq.14) then
         ss2 = 26.4D0*(s**(-.34D0))+59.3D0*(s**.095D0)
        endif
        cs_multidiff = (1.-exp(-ss1))*ss2
        cs_multi = 0.89D0*cs_multidiff

c****************************
c DIFFRACTIVE SCATTERING:
c****************************  

        ss1 = ((x-.85D0)**.75D0)/.64D0
        ss2 = 74.1D0*(x**(-.44D0))+62.D0*(s**.08D0)
        cs_tmp = 0.96D0*(1.D0-exp(-ss1))*ss2
        cross_diffr1 = 0.14D0*cs_tmp
        cross_diffr2 = 0.013D0*cs_tmp
        cs_delta = cross_frag2 - (cross_diffr1+cross_diffr2-cross_diffr)
        if(cs_delta.lt.0.D0) then
          cross_frag2 = 0.D0
          cs_multi = cs_multi+cs_delta
        else
          cross_frag2 = cs_delta
        endif
        cross_diffr = cross_diffr1 + cross_diffr2
        cs_multidiff = cs_multi + cross_diffr

       else
        cross_diffr = 0.D0
        cross_diffr1 = 0.D0
        cross_diffr2 = 0.D0
        cs_multidiff = 0.D0
        cs_multi = 0.D0
       endif

       if (NDIR.eq.3) then

        crossection = cross_res+cross_dir+cs_multidiff+cross_frag2
        RETURN

       else if (NDIR.eq.0) then

        crossection = cross_res+cross_dir+cross_diffr+cross_frag2
        RETURN

       else if (NDIR.eq.2) then

        crossection = cross_res+cross_dir
        RETURN

       else if (NDIR.eq.1) then

        crossection = cross_res
        RETURN

       else if (NDIR.eq.4) then

        crossection = cross_dir
        RETURN

       else if (NDIR.eq.5) then

        crossection = cs_multi
        RETURN

       else if (NDIR.eq.6) then

        crossection = cross_res+cross_dir2
        RETURN

       else if (NDIR.eq.7) then

        crossection = cross_res+cross_dir1
        RETURN

       else if (NDIR.eq.8) then

        crossection = cross_res+cross_dir+cross_diffr1
        RETURN

       else if (NDIR.eq.9) then

        crossection = cross_res+cross_dir+cross_diffr
        RETURN

       else if (NDIR.eq.10) then

        crossection = cross_diffr
        RETURN

       else if ((NDIR.ge.11).and.(NDIR.le.19)) then

        crossection = sig_res(NDIR-10)
        RETURN

       else

        print*,'wrong input NDIR in crossection.f !'
        STOP

       endif
      
       END
