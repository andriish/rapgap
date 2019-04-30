

      SUBROUTINE check_event(Ic,Esum,PXsum,PYsum,PZsum,IQchr,IQbar,Irej)
C***********************************************************************
C
C  check quantum number conservation of final state
C
C                                                (R.E. 08/98)
C
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      COMMON /S_RUN/ SQS, S, Q2MIN, XMIN, ZMIN, kb, kt, a1, a2, Nproc
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
       COMMON /S_CSYDEC/ CBR(102), IDB(49), KDEC(612), LBARP(49)
      COMMON /S_CHP/ S_LIFE(49), ICHP(49), ISTR(49), IBAR(49)
      COMMON /S_MASS1/ AM(49), AM2(49)
      COMMON /S_CNAM/ NAMP (0:49)
      CHARACTER*6 NAMP
      SAVE

      px = 0.D0
      py = 0.D0
      pz = 0.D0
      ee = 0.D0
      ichar = 0
      ibary = 0
      Iprint = 0
      
      PLscale = Esum
      PTscale = 1.D0

      do j=1,np
        l1 = abs(LLIST(J))
        l = mod(llist(j),10000)
        if(l1.lt.10000) then
          px = px + P(j,1)
          py = py + P(j,2)
          pz = pz + P(j,3)
          ee = ee + P(j,4)
          write(6,*)' j,P = ',j,(P(j,i),i=1,4)
c          xm2 = (P(j,4)-P(j,3))*(P(j,4)+P(j,3))
c     &         -P(j,1)**2-P(j,2)**2
c          if(ABS(xm2-P(j,5)**2)/MAX(AM(l1),1.D0).gt.0.1D0) then
c            print *,' vector has incorrect mass (line, particle)',J,l1
c            print *,' mass (theory, actual, p(l,5)): ',
c     &        AM(l1),sqrt(xm2),P(j,5)
c            Iprint = 1
c          endif
          ichar = ichar+sign(1,l)*ICHP(iabs(l))
          ibary = ibary+sign(1,l)*IBAR(iabs(l))
        endif
      enddo

      if(ichar.ne.IQchr) then
        print *,' charge conservation violated',Ic,ichar,iqchr
        Iprint = 1
      endif
      if(ibary.ne.IQbar) then
        print *,' baryon number conservation violated',Ic,ibary,IQbar
        Iprint = 1
      endif
      if(abs((px-PXsum)/MAX(PXsum,PTscale)).gt.1.D-3) then
        print *,' x momentum conservation violated',Ic,px,PXsum
        Iprint = 1
      endif
      if(abs((py-PYsum)/MAX(PYsum,PTscale)).gt.1.D-3) then
        print *,' y momentum conservation violated',Ic,py,PYsum
        Iprint = 1
      endif
      if(abs((pz-Pzsum)/MAX(ABS(PZsum),PLscale)).gt.1.D-3) then
        print *,' z momentum conservation violated',Ic,pz,Pzsum
        Iprint = 1
      endif
      if(abs((ee-Esum)/MAX(Esum,1.D0)).gt.1.D-3) then
        print *,' energy conservation violated',Ic,ee,Esum
        Iprint = 1
      endif

      if(Iprint.ne.0) call print_event(1)

      Irej = Iprint

      END
