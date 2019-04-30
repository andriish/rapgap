

      SUBROUTINE print_event(Iout)
C*********************************************************************
C
C  print final state particles
C
C                                                  (R.E. 03/98)
C
C**********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      COMMON /S_RUN/ SQS, S, Q2MIN, XMIN, ZMIN, kb, kt, a1, a2, Nproc
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
       COMMON /S_CSYDEC/ CBR(102), IDB(49), KDEC(612), LBARP(49)
      COMMON /S_CHP/ S_LIFE(49), ICHP(49), ISTR(49), IBAR(49)
      COMMON /S_MASS1/ AM(49), AM2(49)
      COMMON /S_CNAM/ NAMP (0:49)
      CHARACTER*6 NAMP
      CHARACTER CODE*18
      SAVE

      if(iout.gt.0) then
       
        print *,' --------------------------------------------------'

        if(Nproc.eq.1) then
           print *,' diffractive rho-0 production',Nproc
        else if(Nproc.eq.2) then
           print *,' direct interaction 1',Nproc
        else if(Nproc.eq.3) then
           print *,' direct interaction 2',Nproc
        else if(Nproc.eq.4) then
           print *,' diffractive omega production',Nproc
        else if(Nproc.eq.0) then
           print *,' multi-pion/fragmentation contribution',Nproc
        else if((Nproc.gt.10).and.(Nproc.lt.20)) then
           print *,' resonance production and decay',Nproc-10
        else
           print *,' unknown process',Nproc
        endif

        i0 = 0
        px = 0.D0
        py = 0.D0
        pz = 0.D0
        ee = 0.D0
        ichar = 0
        ibary = 0
        do j=1,np
          l1 = abs(LLIST(J))
          l = mod(llist(j),10000)
          if(l1.lt.10000) then
            px = px + P(j,1)
            py = py + P(j,2)
            pz = pz + P(j,3)
            ee = ee + P(j,4)
            ichar = ichar+sign(1,l)*ICHP(iabs(l))
            ibary = ibary+sign(1,l)*IBAR(iabs(l))
          endif
          if((l1.lt.10000).or.(Iout.GE.2)) then
            i0 = i0+1
            code = '                  '
            code(1:6) = namp(iabs(l))
            if (l .lt. 0) code(7:9) = 'bar'
            write (6,120) i0,CODE,l1*sign(1,l),sign(1,l)*ICHP(iabs(l)),
     &        (P(j,k),k=1,4)
          endif
        enddo
        write (6,122) '   sum: ',px,py,pz,ee
        print *,' charge QN: ',ichar,'    baryon QN: ',ibary
        print *,' --------------------------------------------------'
120     FORMAT(1X,I4,1X,A18,1X,I6,1X,I2,1X,2(F9.3,2X),2(E9.3,2X))
122     FORMAT(7X,A8,20X,2(F9.3,2X),2(E9.3,2X))

      endif

      END
