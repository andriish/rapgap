      SUBROUTINE RGUPEVNT
 
      Implicit None
 
#include "rglujets.inc"
#include "rgpara.inc"
#include "rgpara1.inc"
#include "rgluco.inc"
#include "rgherup.inc"
C...RAPGAP process event common block.
c
c lines marked with special are needed for HERWIG running....
c
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/

      Integer MINT
      Double Precision VINT
      COMMON/PYINT1/MINT(400),VINT(400)

      Integer kb1,kb2,kb3,kdq1,kdq2,kdq3,kd1,kd2,kd3,kq
      Integer MINT111
      Integer nupp
	Integer I,J,ncall,iorg,idiq
      Double Precision psum(4),boost(3),Mtot,Dot
	Data ncall/0/
	Integer Icol
      integer k3(100)
      External Dot
      
	ncall = ncall+1
ccc      write(6,*) ' starting rgupevnt ', ncall,nf1,nf2,nia1,nia2,n
ccc      call pylist(1)
      MINT111=MINT(111)
      If(k(nf1,1).ne.1) k(nf1,1)=-1
      If(k(nf2,1).ne.1) k(nf2,1)=-1
      if(Iabs(K(2,2)).gt.1000) then
      k(nia2,1)=-1
      k(nia2+1,1)=-1
      endif
      if(iabs(K(1,2)).gt.1000) then
      k(nia1,1)=-1
      k(nia1+1,1)=-1
      endif
      call pyedit(11)
      
c      write(6,*) ' after pyedit ', ncall,nf1,nf2,nia1,nia2,n
c      call pylist(2)
c checking diquark and replace it with incoming quark
c for finding it back, incoming quark is 5000+ID      
      kq=-9999
      idiq = -9999
      Do I=3,N
cc        if(iabs(K(I,2)).gt.10000) then
c special treatment when already hadron clustered 
cc           kdq3=1
cc           goto 345
cc        endif
c        if(iabs(K(I,2)).gt.1000.and.iabs(K(I,2)).lt.5000) then
        if(iabs(K(I,2)).gt.1000) then
c          write(6,*) ' diquark found ... ', K(I,2),I,K(I,3)
          iorg = K(I,3)
          KB1=MOD(K(iorg,2)/1000,10)
          KB2=MOD(K(iorg,2)/100,10)
          KB3=MOD(K(iorg,2)/10,10)
          
          KDQ1=MOD(K(I,2)/1000,10)
          KDQ2=MOD(K(I,2)/100,10)
          KDQ3=MOD(K(I,2)/10,10)
          kd1=kb1-kdq1
          kd2=kb2-kdq2
          kd3=kb3-kdq3
c          write(6,*) ' kdq3 ',kdq3,kb3,K(I,2)
          if(kdq3.ne.0) goto 345
c          write(6,*) ' origin ',iorg,K(iorg,2)
c          write(6,*) ' kb1 ',kb1,kb2,kb3
c          write(6,*) ' kdq1 ',kdq1,kdq2,kdq3
c          write(6,*) ' kd1 ',kd1,kd2,kd3
c          if(kd1.ne.0) write(6,*) 'rgupevnt: kd1 ',kd1 
          if(kd2.eq.0.and.kd3.ne.0) then 
            kq = kd3
          elseif(kd2.ne.0.and.kd3.ne.0) then 
c            write(6,*) 'rgupevnt: kd2 ',kd2 
            kq = kd2+kd3
          elseif(kd2.ne.0.and.kd3.eq.0) then 
c            write(6,*) 'rgupevnt: kd2 ',kd2 
            kq = kd2
          endif
c          write(6,*) ' kq = ',kq
          if(kq.eq.-9999)  write(6,*) ' ERROR kq = ',kq
          idiq = i
          Do J=1,3
          P(I,J) = P(iorg,J) - P(I,J)
          Enddo
cold          P(I,5) = PYMASS(kq)
          P(I,5) = 0
          P(I,4) = Dsqrt(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2)
          K(I,2) = kq+isign(1,K(iorg,2))*5000
c          write(6,*) 'kq ... ',kq,K(I,2)
        endif        
      End do
c      write(6,*) ' incoming part 1st method ',(P(idiq,i),i=1,4)
      do j=1,4
        psum(j)=0d0
      end do
      Do I=1,N
      if(I.ne.idiq) then
        if(K(I,1).eq.1.or.K(I,1).eq.2) then
          do j=1,4
            psum(j) = psum(j) + P(i,j)
          end do
        endif
      endif
      End do
      do i =1,3
      boost(i) = psum(i)/psum(4)
      End do
      Mtot = dot(psum,psum)
      CALL pyrobo(0,0,0.D0,0.D0,-boost(1),-boost(2),-boost(3))
c      SPHI(2)=pyangl(P(NIPH,1),P(NIPH,2))
c      CALL pyrobo(0,0,0.D0,-SPHI(2),0.D0,0.D0,0.D0)
c      STHETA(2)=pyangl(P(NIPH,3),P(NIPH,1))
c      CALL pyrobo(0,0,-STHETA(2),0.D0,0.D0,0.D0,0.D0)
c      write(6,*) ' after boost '
c      call pylist(1)
cnew      P(idiq,4) = (Mtot+P(idiq,5)**2)/2d0/dsqrt(Mtot)
cnew      P(idiq,3) = dsqrt(P(idiq,4)**2-P(idiq,5)**2)
cnew      P(idiq,1) = 0
cnew      P(idiq,2) = 0
cnew      P(idiq,5) = 0
      CALL pyrobo(0,0,0.D0,0.D0,boost(1),boost(2),boost(3))
c      write(6,*) ' incoming part 2nd method ',(P(idiq,i),i=1,4)
      
      if(P(idiq,4)/P(2,4).gt.0.7) then
c         write(6,*) ' large x event ',P(idiq,4)/P(2,4)
c         write(6,*) ' Mtot ',dsqrt(Mtot)
c         return
      endif 
      
      NUP=0
	NUP=NUP+1
	ISTUP(NUP) = -1
	IDUP(NUP) = K(2,2)
	ICOLUP(1,NUP) = 0
	ICOLUP(2,NUP) = 0
	Do J = 1,5
	 PUP(J,NUP) = P(2,J)
	Enddo
	NUP=NUP+1
	ISTUP(NUP) = -1
	IDUP(NUP) = K(1,2)
	ICOLUP(1,NUP) = 0
	ICOLUP(2,NUP) = 0
	Do J = 1,5
	 PUP(J,NUP) = P(1,J)
	Enddo
	Icol = 0
      Do I=3,N
	 IF(K(I,1).eq.1.or.K(I,1).eq.2) then
c       write(6,*) ' start loop ',i,k(i,2),N
	    NUP=NUP+1
	    ICOLUP(1,NUP) = 0
	    ICOLUP(2,NUP) = 0
          MOTHUP(1,nup)=3	
          MOTHUP(2,nup)=1	
	    ISTUP(NUP) = 1
	    IDUP(NUP) = K(I,2)
	    if(K(i,1).eq.2) then 
	      Icol = Icol+1
		if(K(i,2).gt.0.and.iabs(K(i,2)).lt.10) then
 	         ICOLUP(1,NUP) = 500+icol
		   ICOLUP(2,NUP) = 0
		   elseif(K(i,2).gt.0.and.iabs(K(i,2)).gt.1000) then
		   ICOLUP(1,NUP) = 500+icol
  	         ICOLUP(2,NUP) = 0
c               write(6,*)  ' parton A K',K(I,2)
               K(I,2) = K(I,2)-isign(1,K(I,2))*5000
	         IDUP(NUP) = K(I,2)
     	         ISTUP(NUP) = -1 
               if(K(I,3).eq.1) then
                 MOTHUP(1,nup) = 2
                 elseif(K(I,3).eq.2) then
                 MOTHUP(1,nup) = 1
                 endif
c               write(6,*) ' parton  A ',IDUP(NUP),K(i,2)
		   elseif(K(i,2).lt.0.and.iabs(K(i,2)).lt.10) then
  	         ICOLUP(1,NUP) = 0
		   ICOLUP(2,NUP) = 500+icol
		   elseif(K(i,2).lt.0.and.iabs(K(i,2)).gt.1000) then
		   ICOLUP(1,NUP) = 0
 	         ICOLUP(2,NUP) = 500+icol
c               write(6,*)  ' parton B K',K(I,2)
               K(I,2) = K(I,2)-isign(1,K(I,2))*5000
	         IDUP(NUP) = K(I,2)
	         ISTUP(NUP) = -1
               if(K(I,3).eq.1) then
                 MOTHUP(1,nup) = 2
                 elseif(K(I,3).eq.2) then
                 MOTHUP(1,nup) = 1
                 endif
c               write(6,*) ' parton B ',IDUP(NUP),K(I,2)
		   elseif(K(i,2).eq.21) then
c look for gluons with KF=21
		      if(ICOLUP(1,NUP-1).gt.ICOLUP(2,NUP-1)) then
			ICOLUP(1,NUP) = 500+icol
			ICOLUP(2,NUP) = 500+icol-1
                  if(MOTHUP(1,NUP-1).eq.1.or.MOTHUP(1,NUP-1).eq.2)then
			ICOLUP(2,NUP) = 500+icol
			ICOLUP(1,NUP) = 500+icol-1
                  endif
		      elseif(ICOLUP(1,NUP-1).lt.ICOLUP(2,NUP-1)) then
			ICOLUP(1,NUP) = 500+icol-1
			ICOLUP(2,NUP) = 500+icol
                  if(MOTHUP(1,NUP-1).eq.1.or.MOTHUP(1,NUP-1).eq.2)then
			ICOLUP(2,NUP) = 500+icol-1
			ICOLUP(1,NUP) = 500+icol
                  endif
			endif
		   elseif(K(i,2).gt.100) then
		   else
  	         ICOLUP(1,NUP) = 0
		   ICOLUP(2,NUP) = 0
		endif   
	    elseif(K(i,1).eq.1) then
c            write(6,*) ' here 1 ',i,k(i,2)
c            write(6,*) ' here 1 col ',ICOLUP(1,NUP-1),ICOLUP(2,NUP-1)
            if(iabs(K(i,2)).gt.1000) then
c               write(6,*) ' incoming A K',k(i,2),nup
               K(I,2) = K(I,2)-isign(1,K(I,2))*5000
c               write(6,*) ' incoming 2nd A K',k(i,2),nup
	         IDUP(NUP) = K(I,2)
	         ISTUP(NUP) = -1
               if(K(i,2).gt.0) then 
		   ICOLUP(1,NUP) = 500+icol
  	         ICOLUP(2,NUP) = 0
               else
  	         ICOLUP(1,NUP) = 0
		   ICOLUP(2,NUP) = 500+icol
               endif 
               MOTHUP(1,nup) = K(I,3)
               if(K(I,3).eq.1) then
                 MOTHUP(1,nup) = 2
                 elseif(K(I,3).eq.2) then
                 MOTHUP(1,nup) = 1
                 endif
c               write(6,*) ' this is for incoming part A ',IDUP(NUP)
	      elseif(ICOLUP(1,NUP-1).ne.0.or.ICOLUP(2,NUP-1).ne.0) then
 		    if(ICOLUP(1,NUP-1).gt.ICOLUP(2,NUP-1)) then
			ICOLUP(1,NUP) = 0
			ICOLUP(2,NUP) = 500+icol
		    elseif(ICOLUP(1,NUP-1).lt.ICOLUP(2,NUP-1)) then
			ICOLUP(1,NUP) = 500+icol
			ICOLUP(2,NUP) = 0
	          endif
	      else
            endif
	    endif
       elseif(K(I,2).LT.1000) then
	    NUP=NUP+1
	    ISTUP(NUP) = 3
	    ICOLUP(1,NUP) = 0
	    ICOLUP(2,NUP) = 0
	    IDUP(NUP) = K(I,2)
	    MOTHUP(2,NUP) = 0  ! special
	    if(K(I,3).eq.1) then 
	      MOTHUP(1,NUP) = 2
	    elseif(K(I,3).eq.2) then 
	      MOTHUP(1,NUP) = 1
	    else
	      MOTHUP(1,NUP) = 3
	      MOTHUP(2,NUP) = 1
	    endif
	 Else
          goto 111
       Endif
	 VTIMUP(NUP) = V(I,5)
	 SPINUP(NUP) = 9

	 Do J = 1,5
	   PUP(J,NUP) = P(I,J)
	 Enddo
	
111   Continue
	Enddo 
      
c      write(6,*) 'before reshuffling'
c      call pylist(7)
      do I=3,NUP
c      write(6,*) ' checking here ',ISTUP(I)
      if(ISTUP(I).lt.0) then 
        iorg=mothup(1,i)
	  Do J = 1,5
	   PUP(J,iorg) = PUP(J,I) 
         PUP(J,I) = 0
	  Enddo
	  ISTUP(iorg) = ISTUP(I)
	  IDUP(iorg) =  IDUP(I)
	  IDUP(iorg) =  IDUP(I)
	  ICOLUP(1,iorg) = ICOLUP(1,I) 
	  ICOLUP(2,iorg) = ICOLUP(2,I) 
        
	  ISTUP(I) = 3
	  IDUP(I) =  0
	  ICOLUP(1,I) = 0
	  ICOLUP(2,I) = 0
	  MOTHUP(1,1) =  0
 	  MOTHUP(2,1) =  0
	  MOTHUP(1,2) =  0
 	  MOTHUP(2,2) =  0
      Endif
      enddo
c      write(6,*) 'before idup'
c      call pylist(7)
     
      nupp =0
      do i=1,nup
      if(idup(i).eq.0) then
      else
         nupp=nupp+1
         istup(nupp)=istup(i)
         icolup(1,nupp) = icolup(1,i)	
         icolup(2,nupp) = icolup(2,i)	
	   IDUP(nupp) =  IDUP(I)
	   MOTHUP(1,nupp) = MOTHUP(1,i)  
 	   MOTHUP(2,nupp) = MOTHUP(2,i)  
	   Do J = 1,5
	     PUP(J,nupp) = PUP(J,I) 
	   Enddo
      endif
      enddo
      nup=nupp
            
345   if(kdq3.ne.0)  then
c here treat event which have a clustered hadron separatly
         MINT(111)=12 
         NUP=0
	   NUP=NUP+1
	   ISTUP(NUP) = -1
	   IDUP(NUP) = K(2,2)
	   ICOLUP(1,NUP) = 0
	   ICOLUP(2,NUP) = 0
	   Do J = 1,5
	    PUP(J,NUP) = P(2,J)
	   Enddo
	   NUP=NUP+1
	   ISTUP(NUP) = -1
	   IDUP(NUP) = K(1,2)
	   ICOLUP(1,NUP) = 0
	   ICOLUP(2,NUP) = 0
	   Do J = 1,5
	     PUP(J,NUP) = P(1,J)
	   Enddo
	   Icol = 0
         Do I=3,N
	   IF(K(I,1).eq.1.or.K(I,1).eq.2) then
	      NUP=NUP+1
	      ICOLUP(1,NUP) = 0
	      ICOLUP(2,NUP) = 0
            MOTHUP(1,nup)=0
            MOTHUP(2,nup)=0	
	      ISTUP(NUP) = 1
	      IDUP(NUP) = K(I,2)
	      if(K(I,3).eq.1) then 
	         MOTHUP(1,NUP) = 2
	      elseif(K(I,3).eq.2) then 
	         MOTHUP(1,NUP) = 1
	      else
	         MOTHUP(1,NUP) = K(I,3)
	      endif
	      MOTHUP(2,NUP) = 0
	    
            MOTHUP(1,nup)=1 ! special
            MOTHUP(2,nup)=1 ! special	
	    
 	      ICOLUP(1,NUP) = 0
	      ICOLUP(2,NUP) = 0
	      if(K(i,1).eq.2) then 
               MOTHUP(1,nup)=1 ! special
               MOTHUP(2,nup)=2 ! special
	         Icol = Icol+1
		   if(K(i,2).gt.0.and.iabs(K(i,2)).lt.10) then
 	            ICOLUP(1,NUP) = 500+icol
		      ICOLUP(2,NUP) = 0
		   elseif(K(i,2).gt.0.and.iabs(K(i,2)).gt.1000) then
  	            ICOLUP(1,NUP) = 0
		      ICOLUP(2,NUP) = 500+icol
		   elseif(K(i,2).lt.0.and.iabs(K(i,2)).lt.10) then
  	            ICOLUP(1,NUP) = 0
		      ICOLUP(2,NUP) = 500+icol
		   elseif(K(i,2).lt.0.and.iabs(K(i,2)).gt.1000) then
 	            ICOLUP(1,NUP) = 500+icol
		      ICOLUP(2,NUP) = 0
		   elseif(K(i,2).eq.21) then
		      if(ICOLUP(1,NUP-1).gt.ICOLUP(2,NUP-1)) then
			   ICOLUP(1,NUP) = 500+icol
			   ICOLUP(2,NUP) = 500+icol-1
		      elseif(ICOLUP(1,NUP-1).lt.ICOLUP(2,NUP-1)) then
			   ICOLUP(1,NUP) = 500+icol-1
			   ICOLUP(2,NUP) = 500+icol
			endif
		   elseif(K(i,2).gt.100) then
		   else
  	            ICOLUP(1,NUP) = 0
		      ICOLUP(2,NUP) = 0
		   endif   
	       elseif(K(i,1).eq.1) then
	         if(ICOLUP(1,NUP-1).ne.0.or.ICOLUP(2,NUP-1).ne.0) then
 		      if(ICOLUP(1,NUP-1).gt.ICOLUP(2,NUP-1)) then
			   ICOLUP(1,NUP) = 0
			   ICOLUP(2,NUP) = 500+icol
		      elseif(ICOLUP(1,NUP-1).lt.ICOLUP(2,NUP-1)) then
			   ICOLUP(1,NUP) = 500+icol
			   ICOLUP(2,NUP) = 0
	            endif
                  MOTHUP(1,nup)=1
                  MOTHUP(2,nup)=2	
	         endif
	       endif
	    else
	       NUP=NUP+1
	       ISTUP(NUP) = 3
	       ICOLUP(1,NUP) = 0
	       ICOLUP(2,NUP) = 0
	       IDUP(NUP) = K(I,2)
	       if(K(I,3).eq.1) then 
	          MOTHUP(1,NUP) = 2
	       elseif(K(I,3).eq.2) then 
	          MOTHUP(1,NUP) = 1
	       else
	          MOTHUP(1,NUP) = K(I,3)
	       endif
	       MOTHUP(2,NUP) = 0  ! special
             MOTHUP(1,nup)=1 ! special
             MOTHUP(2,nup)=1 ! special	
	    Endif
	    VTIMUP(NUP) = V(I,5)
	    SPINUP(NUP) = 9

     	    Do J = 1,5
     	       PUP(J,NUP) = P(I,J)
     	    Enddo
	
	  Enddo 
      endif
      
	IDPRUP = IPRO
      XWGTUP = 1
	SCALUP = 0.1
	AQEDUP = -1.
	AQCDUP = -1.
      
      
c	if(ncall.lt.5) then 
c      write(6,*) ' before print'
c      call pylist(7)
c      endif
      IF(ILHA.ge.10) then        
c         write(6,*) ' starting pyevnt ',ncall,nup
c         write(6,*) ' check VINT(1),VINT(2) ',VINT(1),VINT(2)
c         write(6,*) ' check VINT(141),VINT(142) ',VINT(141),VINT(142)
c         write(6,*) '  x event ',P(idiq,4)/P(2,4)
c	   if(ncall.lt.150) call pylist(7)
	   do i=1,100
	   do j=1,5
	   p(i,j)=0
	   k(i,j)=0
	   enddo
	   enddo
	   CALL PYEVNT
c         write(6,*) ' after pyevnt'
	   if(ncall.lt.5) call pylist(1)
	else
c      write(6,*) ' cascade RGupevnt '
	   Write(21) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &   (IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I),ICOLUP(1,I),
     &   ICOLUP(2,I),(PUP(J,I),J=1,5),VTIMUP(I),SPINUP(I),I=1,NUP)
	
c      write(6,*) ' cascade RGupevnt ',ncall
	   if(ncall.lt.5) call pylist(7)
c	   if(ncall.lt.5) call pylist(1)
	endif
c	call pylist(1)
      MINT(111)=MINT111
      RETURN
      END 
 
