

      SUBROUTINE gamma_h(Ecm,ip1,Imode,ifbad)
C**********************************************************************
C
C     simple simulation of low-energy interactions (R.E. 03/98)
C
C     changed to simulate superposition of reggeon and pomeron exchange 
C     interface to Lund / JETSET 7.4 fragmentation
C                                                  (R.E. 08/98)
C
C     
C
C     input: ip1    incoming particle
C                   13 - p
C                   14 - n
C            Ecm    CM energy in GeV
C            Imode  interaction mode
C                   0 - multi-pion fragmentation
C                   5 - fragmentation in resonance region
C                   1 - quasi-elastic / diffractive interaction 
C                       (p/n-gamma  --> n/p rho)
C                   4 - quasi-elastic / diffractive interaction 
C                       (p/n-gamma  --> n/p omega)
C                   2 - direct interaction (p/n-gamma  --> n/p pi)
C                   3 - direct interaction (p/n-gamma  --> delta pi)
C            IFBAD control flag
C                  (0  all OK,
C                   1  generation of interaction not possible)
C
C**********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      COMMON /S_RUN/ SQS, S, Q2MIN, XMIN, ZMIN, kb, kt, a1, a2, Nproc
      COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
      COMMON /S_CHP/ S_LIFE(49), ICHP(49), ISTR(49), IBAR(49)
      COMMON /S_MASS1/ AM(49), AM2(49)
      COMMON /S_CFLAFR/ PAR(8)
      SAVE

      DIMENSION P_dec(10,5), P_in(5)
      DIMENSION xs1(2), xs2(2), xmi(2), xma(2)
      DIMENSION LL(10), Ijoin(4)

      DOUBLE PRECISION PA1(4), PA2(4), P1(4), P2(4)

      DATA Ic / 0 /

C  second particle is always photon
      IP2 = 1
C  parameters of pi0 suppression
C... 50% suppression:
      a1 = 0.5D0
      a2 = 0.5D0
C... no suppression:
c      a1 = 0.D0
c      a2 = 1.D0
C  parameter of strangeness suppression
      PAR(2) = 0.18D0

      ifbad = 0
      SQS = Ecm
      S = SQS*SQS
      Ic = Ic+1

*     xm_p = AM(13)
*     xm_pi = 0.135
*     plab = (Ecm**2-xm_p**2)/(2.*xm_p)


      IF((Imode.eq.1).or.(Imode.eq.4)) THEN

C***********************************************************************

C  simulation of diffraction

        ipa = ip1
        ipb = ip2

        if(Imode.eq.1) then
          Nproc = 1
          if(ip1.eq.1) then
            ipa = 27
          else if(ip2.eq.1) then
            ipb = 27
          endif
        else if(Imode.eq.4) then
          Nproc = 4
          if(ip1.eq.1) then
            ipa = 32
          else if(ip2.eq.1) then
            ipb = 32
          endif
        endif

        am_a = AM(ipa)
        am_b = AM(ipb)
        if(am_a+am_b.ge.Ecm-1.D-2) am_a = Ecm - am_b-1.D-2

C  find t range
        e1 = 0.5D0*(Ecm + AM(ip1)**2/Ecm - AM(ip2)**2/Ecm)
        if(e1.gt.100.D0*AM(ip1)) then
          pcm1 = e1 - 0.5D0*AM(ip1)**2/e1
        else
          pcm1 = sqrt((e1-AM(ip1))*(e1+AM(ip1)))
        endif
        e3 = 0.5D0*(Ecm + am_a**2/Ecm - am_b**2/Ecm)
        if(e3.gt.100.D0*am_a) then
          pcm3 = e3 - 0.5D0*am_a**2/e3
        else
          pcm3 = sqrt((e3-am_a)*(e3+am_a))
        endif
        t0 = ((AM(ip1)**2-am_a**2-AM(ip2)**2+am_b**2)/(2.D0*Ecm))**2
     &      -(pcm1-pcm3)**2-0.0001D0
        t1 = ((AM(ip1)**2-am_a**2-AM(ip2)**2+am_b**2)/(2.D0*Ecm))**2
     &      -(pcm1+pcm3)**2+0.0001D0

C  sample t
        b = 6.D0
        t = 1.D0/b*log((exp(b*t0)-exp(b*t1))*RNDM(0)+exp(b*t1))

C  kinematics
        pl = (2.D0*e1*e3+t-AM(ip1)**2-am_a**2)/(2.D0*pcm1)
        pt = (pcm3-pl)*(pcm3+pl)
        if(pt.lt.0.D0) then
          pl = sign(pcm3,pl)
          pt = 1.D-6
        else
          pt = sqrt(pt)
        endif
        phi = 6.28318530717959D0*RNDM(0)

        LLIST(1) = ipa
        P(1,4) = e3
        P(1,1) = SIN(phi)*pt
        P(1,2) = COS(phi)*pt
        P(1,3) = pl
        P(1,5) = am_a
        LLIST(2) = ipb
        P(2,1) = -P(1,1)
        P(2,2) = -P(1,2)
        P(2,3) = -P(1,3)
        P(2,4) = Ecm - P(1,4)
        P(2,5) = am_b
        np = 2

        call DECSIB

      ELSE IF((Imode.eq.2).or.(Imode.eq.3)) THEN

C***********************************************************************

C  simulation of direct p-gamma process

        if(ip1.eq.13) then
C  projectile is a proton
          if(Imode.eq.2) then
            Nproc = 2
            ipa = 14
            ipb = 7
          else
            Nproc = 3
            if(rndm(0).gt.0.25) then
              ipa = 40
              ipb = 8
            else
              ipa = 42
              ipb = 7
            endif
          endif
        else if(ip1.eq.14) then
C  projectile is a neutron
          if(Imode.eq.2) then
            Nproc = 2
            ipa = 13
            ipb = 8
          else
            Nproc = 3
            if(rndm(0).gt.0.25) then
              ipa = 43
              ipb = 7
            else
              ipa = 41
              ipb = 8
            endif
          endif
        endif

        am_a = AM(ipa)
        am_b = AM(ipb)
        if(am_a+am_b.ge.Ecm-1.e-3) am_a = Ecm - am_b-1.D-3

C  find t range
        e1 = 0.5D0*(Ecm + AM(ip1)**2/Ecm - AM(ip2)**2/Ecm)
        if(e1.gt.100.D0*AM(ip1)) then
          pcm1 = e1 - 0.5D0*AM(ip1)**2/e1
        else
          pcm1 = sqrt((e1-AM(ip1))*(e1+AM(ip1)))
        endif
        e3 = 0.5D0*(Ecm + am_a**2/Ecm - am_b**2/Ecm)
        if(e3.gt.100.D0*am_a) then
          pcm3 = e3 - 0.5D0*am_a**2/e3
        else
          pcm3 = sqrt((e3-am_a)*(e3+am_a))
        endif
        t0 = ((AM(ip1)**2-am_a**2-AM(ip2)**2+am_b**2)/(2.D0*Ecm))**2
     &      -(pcm1-pcm3)**2-0.0001D0
        t1 = ((AM(ip1)**2-am_a**2-AM(ip2)**2+am_b**2)/(2.D0*Ecm))**2
     &      -(pcm1+pcm3)**2+0.0001D0

C  sample t
        b = 12.D0
        t = 1./b*log((exp(b*t0)-exp(b*t1))*RNDM(0)+exp(b*t1))

C  kinematics
        pl = (2.D0*e1*e3+t-AM(ip1)**2-am_a**2)/(2.D0*pcm1)
        pt = (pcm3-pl)*(pcm3+pl)
        if(pt.lt.0.D0) then
          pl = sign(pcm3,pl)
          pt = 1.D-6
        else
          pt = sqrt(pt)
        endif
        phi = 6.28318530717959D0*RNDM(0)

        LLIST(1) = ipa
        P(1,4) = e3
        P(1,1) = SIN(phi)*pt
        P(1,2) = COS(phi)*pt
        P(1,3) = pl
        P(1,5) = am_a
        LLIST(2) = ipb
        P(2,1) = -P(1,1)
        P(2,2) = -P(1,2)
        P(2,3) = -P(1,3)
        P(2,4) = Ecm - P(1,4)
        P(2,5) = am_b
        np = 2

        call DECSIB

      ELSE

C***********************************************************************

C  simulation of multiparticle production via fragmentation

          Nproc = 0

          SIG_reg  = 129.D0*(S-AM(13)**2)**(-0.4525D0)
          SIG_pom  = 67.7D0*(S-AM(13)**2)**0.0808D0

          if(S.gt.2.6D0) then
            prob_reg = SIG_reg/(SIG_pom+SIG_reg)
          else
            prob_reg = 1.D0
          endif

          ptu =.36D0+.08D0*log10(sqs/30.D0)

          s1 = 1.2D0
          s2 = 0.6D0
          as1 = s1**2/S
          as2 = s2**2/S
          if(s1+s2.ge.sqs-0.2) then
            prob_reg = 1.D0
          endif

          itry = 0
 100      continue
          Istring = 0

C  avoid infinite looping
          itry = itry+1
          if(itry.gt.50) then
            print *,' gamma_h: more than 50 internal rejections,'
            print *,' called with ip1,ip2,Ecm,Imode:',ip1,ip2,Ecm,Imode
            ifbad = 1
            return
          endif

C  simulate reggeon (one-string topology)

          if(RNDM(0).lt.prob_reg) then

            do i=1,1000
              call valences(IP1,Ifl1a,Ifl1b)
              call valences(IP2,Ifl2a,Ifl2b)
              if(Ifl1b.eq.-Ifl2b) goto 200
            enddo
            print *,'gamma_h: simulation of reggeon impossible:',ip1,ip2
            goto 100
            
 200        continue

            np = 0
            Istring = 1

            ee = Ecm/2.D0
 250        continue
              pt = ptu*sqrt(-log(max(1.D-10,RNDM(0))))
            if(pt.ge.ee) goto 250
            phi = 6.2831853D0*RNDM(0)
            px = pt*COS(phi)
            py = pt*SIN(phi)
            
            pz = SQRT(ee**2-px**2-py**2)
            call lund_put(1,Ifl1a,px,py,pz,ee)
            px = -px
            py = -py
            pz = -pz
            call lund_put(2,Ifl2a,px,py,pz,ee)
            Ijoin(1) = 1
            Ijoin(2) = 2
            call pyjoin(2,Ijoin)

            call lund_frag(Ecm,NP)
            if(NP.lt.0) then
              if(Ideb.ge.5) 
     &          print *,' gamma_h: rejection (1) by lund_frag, sqs:',Ecm
              NP = 0
              goto 100
            endif

            do i=1,NP
              call lund_get(i,LLIST(i),
     &                      P(i,1),P(i,2),P(i,3),P(i,4),P(i,5))
            enddo
              

C  simulate pomeron (two-string topology)

          else

            call valences(IP1,Ifl1a,Ifl1b)
            call valences(IP2,Ifl2a,Ifl2b)
            if(Ifl1a*Ifl2a.lt.0) then
              j = Ifl2a
              Ifl2a = Ifl2b
              Ifl2b = j
            endif

            pl1 = (1.D0+as1-as2)
            ps1 = 0.25D0*pl1**2-as1
            if(ps1.le.0.D0) then
              if(Ideb.ge.5) print *,' rejection by x-limits (1) ',Ecm
              prob_reg = 1.D0
              goto 100
            endif
            ps1 = sqrt(ps1)
            xmi(1) = 0.5D0*pl1-ps1
            xma(1) = 0.5D0*pl1+ps1

            pl2 = (1.D0+as2-as1)
            ps2 = 0.25D0*pl2**2-as2
            if(ps2.le.0.D0) then
              if(Ideb.ge.5) print *,' rejection by x-limits (2) ',Ecm
              prob_reg = 1.D0
              goto 100
            endif
            ps2 = sqrt(ps2)
            xmi(2) = 0.5D0*pl2-ps2
            xma(2) = 0.5D0*pl2+ps2

            if((xmi(1).ge.xma(1)+0.05D0).or.
     &         (xmi(2).ge.xma(2)+0.05D0)) then
              if(Ideb.ge.5) print *,' rejection by x-limits (3) ',Ecm
              prob_reg = 1.D0
              goto 100
            endif
            call PO_SELSX2(xs1,xs2,xmi,xma,as1,as2,Irej)
            if(Irej.ne.0) then
              if(Ideb.ge.5) print *,
     &          'gamma_h: rejection by PO_SELSX2, sqs,m1,m2:',Ecm,s1,s2
              prob_reg = 1.D0
              goto 100
            endif

            NP = 0
            Istring = 1

            ee = SQRT(XS1(1)*XS2(1))*Ecm/2.D0
 260        continue
              pt = ptu*sqrt(-log(max(1.D-10,RNDM(0))))
            if(pt.ge.ee) goto 260
            phi = 6.2831853D0*RNDM(0)
            px = pt*COS(phi)
            py = pt*SIN(phi)

            PA1(1) = px
            PA1(2) = py
            PA1(3) = XS1(1)*Ecm/2.D0
            PA1(4) = PA1(3)

            PA2(1) = -px
            PA2(2) = -py
            PA2(3) = -XS2(1)*Ecm/2.D0
            PA2(4) = -PA2(3)

            XM1 = 0.D0
            XM2 = 0.D0
            call PO_MSHELL(PA1,PA2,XM1,XM2,P1,P2)
            px = P1(1)
            py = P1(2)
            pz = P1(3)
            ee = P1(4)
            call lund_put(1,Ifl1a,px,py,pz,ee)
            px = P2(1)
            py = P2(2)
            pz = P2(3)
            ee = P2(4)
            call lund_put(2,Ifl2a,px,py,pz,ee)

            Ijoin(1) = 1
            Ijoin(2) = 2
            call pyjoin(2,Ijoin)

            ee = SQRT(XS1(2)*XS2(2))*Ecm/2.D0
 270        continue
              pt = ptu*sqrt(-log(max(1.D-10,RNDM(0))))
            if(pt.ge.ee) goto 270
            phi = 6.2831853D0*RNDM(0)
            px = pt*COS(phi)
            py = pt*SIN(phi)

            PA1(1) = px
            PA1(2) = py
            PA1(3) = XS1(2)*Ecm/2.D0
            PA1(4) = PA1(3)

            PA2(1) = -px
            PA2(2) = -py
            PA2(3) = -XS2(2)*Ecm/2.D0
            PA2(4) = -PA2(3)

            XM1 = 0.D0
            XM2 = 0.D0
            call PO_MSHELL(PA1,PA2,XM1,XM2,P1,P2)

            px = P1(1)
            py = P1(2)
            pz = P1(3)
            ee = P1(4)
            call lund_put(3,Ifl1b,px,py,pz,ee)
            px = P2(1)
            py = P2(2)
            pz = P2(3)
            ee = P2(4)
            call lund_put(4,Ifl2b,px,py,pz,ee)

            Ijoin(1) = 3
            Ijoin(2) = 4
            call pyjoin(2,Ijoin)

            call lund_frag(Ecm,NP)
            if(NP.lt.0) then
              if(Ideb.ge.5) 
     &          print *,' gamma_h: rejection (2) by lund_frag, sqs:',Ecm
              NP = 0
              prob_reg = prob_reg+0.1D0
              goto 100
            endif

            do i=1,NP
              call lund_get(i,LLIST(i),
     &                      P(i,1),P(i,2),P(i,3),P(i,4),P(i,5))
            enddo
              
          endif

          if(Ideb.ge.10) then
            print *,' multi-pion event',Istring,NP
            call print_event(1)
          endif

c... for fragmentation in resonance region:
          if (Imode.eq.5) goto 400

C  leading baryon/meson effect

          do j=1,np
            if(((LLIST(J).eq.13).or.(LLIST(J).eq.14))
     &         .and.(p(j,3).lt.0.D0)) then
              if(rndm(0).lt.(2.D0*p(j,4)/Ecm)**2) goto 100
            endif
            if((LLIST(J).ge.6).and.(LLIST(J).le.8)
     &         .and.(p(j,3).lt.-0.4D0)) then
              if(rndm(0).lt.(2.D0*p(j,4)/Ecm)**2) goto 100
            endif
          enddo

C  remove elastic/diffractive channels

          ima_0  = 0
          imb_0  = 0
          ima_1  = 0
          imb_1  = 0
          ima_2  = 0
          imb_2  = 0
          imul = 0

          if(ip1.eq.1) then
            iba_0 = 6
            iba_1 = 27
            iba_2 = 32
          else
            iba_0 = ip1
            iba_1 = ip1
            iba_2 = ip1
          endif
          if(ip2.eq.1) then
            ibb_0 = 6
            ibb_1 = 27
            ibb_2 = 32
          else
            ibb_0 = ip2
            ibb_1 = ip2
            ibb_2 = ip2
          endif

          do j=1,np
            l1 = abs(LLIST(J))
            if(l1.lt.10000) then
              if(LLIST(J).eq.iba_0) ima_0 = 1
              if(LLIST(J).eq.iba_1) ima_1 = 1
              if(LLIST(J).eq.iba_2) ima_2 = 1
              if(LLIST(J).eq.ibb_0) imb_0 = 1
              if(LLIST(J).eq.ibb_1) imb_1 = 1
              if(LLIST(J).eq.ibb_2) imb_2 = 1
              imul = imul+1
            endif
          enddo 

          if(imul.eq.2) then
            if(ima_0*imb_0.eq.1) goto 100
            if(ima_1*imb_1.eq.1) goto 100
            if(ima_2*imb_2.eq.1) goto 100
          endif

C  remove direct channels

          if((imul.eq.2).and.
     &       (ip2.eq.1).and.((ip1.eq.13).or.(ip1.eq.14))) then

            ima_0  = 0
            imb_0  = 0
            ima_1  = 0
            imb_1  = 0
            ima_2  = 0
            imb_2  = 0
            ima_3  = 0
            imb_3  = 0

            if(ip1.eq.13) then
              iba_0 = 14
              ibb_0 = 7
              iba_1 = 40
              ibb_1 = 8
              iba_2 = 42
              ibb_2 = 7
              iba_3 = 13
              ibb_3 = 23
            else
              iba_0 = 13
              ibb_0 = 8
              iba_1 = 43
              ibb_1 = 7
              iba_2 = 41
              ibb_2 = 8
              iba_3 = 14
              ibb_3 = 23
            endif
  
            do j=1,np
              l1 = abs(LLIST(J))
              if(l1.lt.10000) then
                if(LLIST(J).eq.iba_0) ima_0 = 1
                if(LLIST(J).eq.iba_1) ima_1 = 1
                if(LLIST(J).eq.iba_2) ima_2 = 1
                if(LLIST(J).eq.iba_3) ima_3 = 1
                if(LLIST(J).eq.ibb_0) imb_0 = 1
                if(LLIST(J).eq.ibb_1) imb_1 = 1
                if(LLIST(J).eq.ibb_2) imb_2 = 1
                if(LLIST(J).eq.ibb_3) imb_3 = 1
              endif
            enddo
            
            if(ima_0*imb_0.eq.1) goto 100
            if(ima_1*imb_1.eq.1) goto 100
            if(ima_2*imb_2.eq.1) goto 100
            if(ima_3*imb_3.eq.1) goto 100

          endif

C  suppress events with many pi0's

          ima_0 = 0
          imb_0 = 0
          do j=1,np
C  neutral mesons
            if(LLIST(J).eq.6) ima_0 = ima_0+1
            if(LLIST(J).eq.11) ima_0 = ima_0+1
            if(LLIST(J).eq.12) ima_0 = ima_0+1
            if(LLIST(J).eq.21) ima_0 = ima_0+1
            if(LLIST(J).eq.22) ima_0 = ima_0+1
            if(LLIST(J).eq.23) ima_0 = ima_0+1
            if(LLIST(J).eq.24) ima_0 = ima_0+1
            if(LLIST(J).eq.27) ima_0 = ima_0+1
            if(LLIST(J).eq.32) ima_0 = ima_0+1
            if(LLIST(J).eq.33) ima_0 = ima_0+1
C  charged mesons
            if(LLIST(J).eq.7) imb_0 = imb_0+1
            if(LLIST(J).eq.8) imb_0 = imb_0+1
            if(LLIST(J).eq.9) imb_0 = imb_0+1
            if(LLIST(J).eq.10) imb_0 = imb_0+1
            if(LLIST(J).eq.25) imb_0 = imb_0+1
            if(LLIST(J).eq.26) imb_0 = imb_0+1
          enddo

          prob_1 = a1*DBLE(imb_0)/max(DBLE(ima_0+imb_0),1.D0)+a2

          if(RNDM(0).GT.prob_1) goto 100


C  correct multiplicity at very low energies

          ND = 0

          E_ref_1 = 1.6D0
          E_ref_2 = 1.95D0

          if((imul.eq.3)
     &       .and.(Ecm.gt.E_ref_1).and.(Ecm.lt.E_ref_2)) then

            ima_0 = 0
            ima_1 = 0
            ima_2 = 0
            imb_0 = 0
            imb_1 = 0
            iba_0 = 0
            iba_1 = 0
            iba_2 = 0
            ibb_0 = 0
            ibb_1 = 0
C  incoming proton
            if(ip1.eq.13) then
              iba_0 = 13
              iba_1 = 7
              iba_2 = 8
              ibb_0 = 14
              ibb_1 = 6
C  incoming neutron
            else if(ip1.eq.14) then
              iba_0 = 14
              iba_1 = 7
              iba_2 = 8
              ibb_0 = 13
              ibb_1 = 6
            endif
            do j=1,np
              if(LLIST(J).eq.iba_0) ima_0 = ima_0+1
              if(LLIST(J).eq.iba_1) ima_1 = ima_1+1
              if(LLIST(J).eq.iba_2) ima_2 = ima_2+1
              if(LLIST(J).eq.ibb_0) imb_0 = imb_0+1
              if(LLIST(J).eq.ibb_1) imb_1 = imb_1+1
            enddo

C  N gamma --> N pi+ pi-
            if(ima_0*ima_1*ima_2.eq.1) then
              Elog = LOG(Ecm)
              Elog_1 = LOG(E_ref_1) 
              Elog_2 = LOG(E_ref_2) 
              prob = 0.1D0*4.D0/(Elog_2-Elog_1)**2
     &                   *(Elog-Elog_1)*(Elog_2-Elog)

              if(RNDM(0).lt.prob) then
                LL(1) = ip1
                LL(2) = 7
                LL(3) = 8
                LL(4) = 6
                ND = 4
              endif

            endif

          endif


          E_ref_1 = 1.95D0
          E_ref_2 = 2.55D0

          if((imul.eq.4)
     &       .and.(Ecm.gt.E_ref_1).and.(Ecm.lt.E_ref_2)) then

            ima_0 = 0
            ima_1 = 0
            ima_2 = 0
            imb_0 = 0
            imb_1 = 0
            iba_0 = 0
            iba_1 = 0
            iba_2 = 0
            ibb_0 = 0
            ibb_1 = 0
C  incoming proton
            if(ip1.eq.13) then
              iba_0 = 13
              iba_1 = 7
              iba_2 = 8
              ibb_0 = 14
              ibb_1 = 6
C  incoming neutron
            else if(ip1.eq.14) then
              iba_0 = 14
              iba_1 = 7
              iba_2 = 8
              ibb_0 = 13
              ibb_1 = 6
            endif
            do j=1,np
              if(LLIST(J).eq.iba_0) ima_0 = ima_0+1
              if(LLIST(J).eq.iba_1) ima_1 = ima_1+1
              if(LLIST(J).eq.iba_2) ima_2 = ima_2+1
              if(LLIST(J).eq.ibb_0) imb_0 = imb_0+1
              if(LLIST(J).eq.ibb_1) imb_1 = imb_1+1
            enddo

C  N gamma --> N pi+ pi- pi0
            if(ima_0*ima_1*ima_2*imb_1.eq.1) then
              Elog = LOG(Ecm)
              Elog_2 = LOG(E_ref_2) 
              Elog_1 = LOG(E_ref_1) 
              prob = 0.1D0*4.D0/(Elog_2-Elog_1)**2
     &                   *(Elog-Elog_1)*(Elog_2-Elog)

              if(RNDM(0).lt.prob) then
                if(ip1.eq.13) then
                  LL(1) = 14
                  LL(2) = 7
                  LL(3) = 7
                  LL(4) = 8
                else
                  LL(1) = 13
                  LL(2) = 7
                  LL(3) = 8
                  LL(4) = 8
                endif
                ND = 4
              endif

            endif

          endif


          if(ND.gt.0) then
            P_in(1) = 0.D0
            P_in(2) = 0.D0
            P_in(3) = 0.D0
            P_in(4) = Ecm
            P_in(5) = Ecm
            call DECPAR(0,P_in,ND,LL,P_dec)
            Iflip = 0
            do j=1,ND
              LLIST(j) = LL(j)
              do k=1,5
                P(j,k) = P_dec(j,k)
              enddo
              if(((LLIST(j).eq.13).or.(LLIST(j).eq.14))
     &           .and.(P(j,3).lt.0.D0)) Iflip = 1
            enddo
            if(Iflip.ne.0) then
              do j=1,ND
                P(j,3) = -P(j,3)
              enddo
            endif
            NP = ND
          endif

c... for fragmentation in resonance region:
  400     continue

          call DECSIB

      ENDIF

      if(Ideb.ge.10) then
        if(Ideb.ge.20) then
          call print_event(2)
        else
          call print_event(1)
        endif
      endif

      IQchr = ICHP(ip1)+ICHP(ip2)
      IQbar = IBAR(ip1)+IBAR(ip2)
      call check_event(-Ic,Ecm,0.D0,0.D0,0.D0,IQchr,IQbar,Irej)

      end
