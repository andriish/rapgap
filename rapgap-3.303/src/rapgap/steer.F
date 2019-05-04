	subroutine steer
	
	Implicit None
	
      character *132 TXT
	character *6 char
	character *10 char2
      character *132 c1
	Integer i1,i2,i3
	Integer inte,i,j
	real r3
	
	
	Integer nloop,nmax
	Parameter (nmax=1000)
      character *4   para
	Real Rval 
	Integer Ival,Ld,Le
      character*132 Cval 
	Common/steering/Nloop,Ld(nmax),Le(nmax),Ival(nmax),
     &	para(nmax),Rval(nmax),Cval(nmax)
	Integer Nevent
	Common/steer1/Nevent
	character *132 filnam
	Integer lunhbk
	Common/steer2/lunhbk,filnam
      LOGICAL         HBKOUT
      COMMON /QHBKLO/ HBKOUT
	
      character *132 hm2filnam
	  Common/steerhm2/hm2filnam
      LOGICAL         hm2OUT
      COMMON /Qhm2LO/ hm2OUT

      character *132 hm3filnam
      integer hm3filfor
	  Common/steerhm3/hm3filnam,hm3filfor
      LOGICAL         hm3OUT
      COMMON /Qhm3LO/ hm3OUT

      character *132 rvtfilnam
      character *132 rvtrun
      character *132 rvtana
	  Integer Nrvtana
	  Common/steerrvt/Nrvtana,rvtfilnam,rvtrun,rvtana(50)
      LOGICAL         rvtOUT,rvtOUT1,rvtOUT2
      COMMON /QrvtLO/ rvtOUT
	
	Nrvtana=0
	nloop = 0
   20 Read(5,101,END=100) TXT
c      write(6,*) ' reading txt: ',txt
      If(TXT(1:1).EQ.'*') then
*          WRITE(6,*)  TXT
         GOTO 20
      Endif
c check for special cards	
      If(TXT(1:6).EQ.'NEVENT') then
         read(txt,*) char,i1
         NEVENT = i1
	   write(6,*) ' Nevents = ',i1
         GOTO 20
      Endif
      If(TXT(1:6).EQ.'HBKOUT') then
         read(txt,*) char,lunhbk,filnam
	   write(6,*) char,'=',lunhbk,' ',filnam
	   HBKOUT=.TRUE.
         GOTO 20
      Endif

      If(TXT(1:6).EQ.'HM2OUT') then
         read(txt,*) char,hm2filnam
	   write(6,*) char,'=',hm2filnam
	   HM2OUT=.TRUE.
         GOTO 20
      Endif
      
      If(TXT(1:6).EQ.'HM3OUT') then
         read(txt,*) char,hm3filnam,hm3filfor
	   write(6,*) char,'=',hm3filnam,' ',hm3filfor
	   HM3OUT=.TRUE.
         GOTO 20
      Endif
      
      If(TXT(1:6).EQ.'RVTOUT') then
       read(txt,*) char,rvtfilnam
       write(6,*) char,'=',rvtfilnam
       write(6,*) char,'=',rvtrun
	   RVTOUT1=.TRUE.
       GOTO 20
      Endif       


    
      If(TXT(1:6).EQ.'RVTANA') then
         if (Nrvtana. LE. 50) then
         Nrvtana=Nrvtana+1
         read(txt,*) char,rvtana(Nrvtana)
	     write(6,*) char,'=',rvtana(Nrvtana)
	     else
         write(*,*)'Too many analyses: ',Nrvtana
         endif
         RVTOUT2=.TRUE.
         GOTO 20
      Endif      
      
      
	If(TXT(1:4).EQ.'END$') goto 90
	inte = 0 
      do j=2,10
	i=j-1
	if(TXT(i:i).EQ."'".and.
     &   (TXT(J:J).EQ.'H'.OR.TXT(J:J).EQ.'I'.OR.TXT(J:J).EQ.'J'.
     &   OR.TXT(J:J).EQ.'K'.OR.TXT(J:J).EQ.'L'.OR.TXT(J:J).EQ.'M'.
     &   OR.TXT(J:J).EQ.'N')) then
         inte = 1
	   elseif(TXT(i:i).EQ."'".and.TXT(J:J).EQ.'U') then
         inte=2
	endif
      if(TXT(i:i).EQ."'".and.TXT(J:J+2).EQ.'PDF') then
         inte=2
	endif
      
	if(TXT(i:i).EQ."'".and.
     &   (TXT(J:J+1).EQ.'KT'.OR.TXT(J:J+1).EQ.'kt')) Then
         inte = 0
      endif
	if(TXT(i:i).EQ."'".and.
     &   (TXT(J:J+2).EQ.'LQ2'.OR.TXT(J:J+2).EQ.'lq2')) Then
         inte = 0
      endif
	enddo
	nloop = nloop + 1
	Ival(nloop) = -9999
	Rval(nloop) = -9999.
	Cval(nloop) = ' '
	if(inte.eq.1) then
c      write(6,*) ' txt ',txt
      read(txt,*) char,i1,i2,i3,char2
c	write(6,*) ' output Integer :',char,i1,i2,i3 
	Ival(nloop) = i3
	elseif(inte.eq.2) then
      read(txt,*) char,c1
c	write(6,*) ' output Charcter :',char,c1,nloop
	Cval(nloop) = c1
      else
      read(txt,*) char,i1,i2,r3,char2
c	write(6,*) ' output real    :',char,i1,i2,r3
	Rval(nloop) = r3
	endif
	Ld(nloop) = i1
	Le(nloop) = i2
	para(nloop) = char
	
c	write(6,*) ' output : ',nloop,para(nloop),Ival(nloop),rval(nloop),cval(nloop)
	
      goto 20
	
   90 Continue
  100 Continue
      RVTOUT=(RVTOUT1 .AND. RVTOUT2)
      Return
	
  101 Format(A132)
      End
