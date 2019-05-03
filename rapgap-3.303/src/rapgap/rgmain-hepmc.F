      PROGRAM RGMAIN
	Implicit None
#include "rgfull.inc"
#include "rgdisdif.inc"
#include "rgluco.inc"
#include "rgludat1.inc"
#include "rgpara.inc"
#include "rgpypars.inc"
#include "rgdiffr.inc"
#include "rgscq2.inc"
#include "rglq2.inc"
#include "rghsunts.inc"
#include "rghsoptn.inc"
#include "rghscuts.inc"
#include "rghsvglp.inc"
#include "rgeffic.inc"
#include "rgpart.inc"
#include "rgstru.inc"
#include "rgrapgki.inc"

      character*512 hepmcout,lhefile
      common/CAHEPMCOUT/hepmcout,lhefile
      integer iev,isub
      double precision xsec,xsece

	Integer Nevent
	Common/steer1/Nevent
	Integer K1,K2,LUX

	real timeleft
      Integer Minuts
      External Minuts
	External pydata
      double precision x1,x2,q2pdfeval,xf1mom,xf2mom,getq2q
      integer pdf1,pdf2

      Integer I,ISEED
      character*255 raseed
      character*255 test
C-----------------------
      character *132 hm2filnam
	  Integer hm2lunhbk
	  Common/steerhm2/hm2lunhbk,hm2filnam
      LOGICAL         hm2OUT
      COMMON /Qhm2LO/ hm2OUT

      character *132 hm3filnam
      character *132 hm3filfor
	  Integer hm3lunhbk
	  Common/steerhm3/hm3lunhbk,hm3filnam,hm3filfor
      LOGICAL         hm3OUT
      COMMON /Qhm3LO/ hm3OUT

      character *132 rvtfilnam
      character *132 rvtana
	  Integer rvtlunhbk
	  Integer rvtanalunhbk
	  Common/steerrvt/rvtlunhbk,rvtfilnam,rvtanalunhbk,rvtana
      LOGICAL         rvtOUT,rvtOUT1,rvtOUT2
      COMMON /QrvtLO/ rvtOUT,rvtOUT1,rvtOUT2
      INTEGER rivetinit,rivetdone
C-----------------------
C...HepMC3
      INTEGER NMXHEP
      PARAMETER (NMXHEP=4000)
      COMMON /HEPEVT/  NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                 JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &                 VHEP(4,NMXHEP)
      INTEGER          NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      DOUBLE PRECISION PHEP,VHEP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C...These variables will be used for conversion. Block size is different
C...in Pythia6 and in HepMC3, so the most simple portable way is to have
C... a second block of same size as in HepMC3 and  copy the content of 
C...block directly.

      INTEGER NMXHEPL,J
      INTEGER new_writer,set_hepevt_address,convert_event,write_event
      INTEGER clear_event,delete_writer
      PARAMETER (NMXHEPL=10000)
      COMMON /HEPEVTL/  NEVHEPL,NHEPL,ISTHEPL(NMXHEPL),IDHEPL(NMXHEPL),
     &           JMOHEPL(2,NMXHEPL),JDAHEPL(2,NMXHEPL),PHEPL(5,NMXHEPL),
     &                 VHEPL(4,NMXHEPL)
      INTEGER          NEVHEPL,NHEPL,ISTHEPL,IDHEPL,JMOHEPL,JDAHEPL
      DOUBLE PRECISION PHEPL,VHEPL
      
      INTEGER OUTID, HEPMC3STATUS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      call dummysub(pydata)
      
C---initialise ARIADNE parameters, now done via block data
C---initialise PYTHIA 6 parameters, via pythia block data 
C     initialize random number generator
c      ISEED = 213123
c      ISEED = Iabs(MINUTS())
c      ISEED = 44788029
      CALL GETENV("RASEED", raseed)
c note RASEED can have only 10 digits
c i.e   export RASEED='1234567890
      read(raseed,'(I30)') iseed
      WRITE (6,*) ' test raseed:',Iseed
      LUX = 4
      K1=0
      K2=0
      CALL RLUXGO(LUX,ISEED,K1,K2)
C---initialise RAPGAP parameters
      CALL GRAINI
C-- read in parameters from file 
      Call Steer
      if (hm3OUT) then
      OUTID=new_writer(0,1,hm3filnam)
      NEVHEP=-123456
      HEPMC3STATUS=set_hepevt_address(NEVHEPL)
      endif
      
      if (rvtOUT1 .AND. rvtOUT2 ) then
      HEPMC3STATUS=rivetinit()
      
      endif
      
      
C-- change standard parameters of RAPGAP 	
	Call rapcha
C-- change standard parameters of HERACLES	
	Call hercha
C-- change standard parameters of JETSET
	Call pytcha
	IF(IFPS.EQ.10) then
C Initialize ARIADNE
         CALL ARINIT('RAPGAP')
C-- change standard parameters of ARIADNE
	   Call aricha
      endif
C--- CALCULATE X SECTION
      CALL PTIME(' rapgap  ',1,0)
      CALL RAPGAP
C--- print x section
      CALL RAEND(1)
C--- event generation
      DO 10 I=1,Nevent
         CALL PTIME(' event  ',1,0)
c	   call timel(timeleft)
c	   if(timeleft.le.10.) then
c	      write(6,*) ' *** time limit reached. time left = ',timeleft
c            goto 123
c	   endif
         CALL EVENT
C--- user analysis routine
c         CALL ANALYS
C--- user analysis routine (hztool)
cc         CALL hzraana
         iev = I
         isub = ipro
c hepmc xsec in pb        
         xsec = avgi *1000.
         xsece = sd * 1000.          

         x1=XEL
         x2=XPR
         q2pdfeval=getq2q();
         xf1mom=XEL*PLEPIN
         xf2mom=XPR*PPIN
         pdf1=MSTP(52)
         pdf2=MSTP(51)
          
         call convhepmc(iev, isub, xsec, xsece,
     &  KE, KP,x1,x2,q2pdfeval,xf1mom,xf2mom ,pdf1,pdf2  )
         CALL PTIME(' event  ',2,0)
         
         if (hm3OUT) then
          NEVHEPL=IEV
          NHEPL=NHEP
           DO 500 J=1,NHEP
          ISTHEPL(J)=ISTHEP(J)
          IDHEPL(J)=IDHEP(J)
          JMOHEPL(1,J)=JMOHEP(1,J)
          JMOHEPL(2,J)=JMOHEP(2,J)
          JDAHEPL(1,J)=JDAHEP(1,J)
          JDAHEPL(2,J)=JDAHEP(2,J)
          PHEPL(1,J)=PHEP(1,J)
          PHEPL(2,J)=PHEP(2,J)
          PHEPL(3,J)=PHEP(3,J)
          PHEPL(4,J)=PHEP(4,J)
          PHEPL(5,J)=PHEP(5,J)
          VHEPL(1,J)=VHEP(1,J)
          VHEPL(2,J)=VHEP(2,J)
          VHEPL(3,J)=VHEP(3,J)
          VHEPL(4,J)=VHEP(4,J)
  500     CONTINUE  
          HEPMC3STATUS=convert_event(OUTID)
C...Note: no XS uncertainty
          HEPMC3STATUS=write_event(OUTID)
          HEPMC3STATUS=clear_event(OUTID)  
         endif
         
C---
   10 CONTINUE
123   Continue
      if (hm3OUT) HEPMC3STATUS=delete_writer(OUTID)
      if (rvtOUT1 .AND. rvtOUT2 ) then
      HEPMC3STATUS=rivetdone()
      
      endif
      
C---PRINT NR OF GENERATED EVENTS
      CALL RAEND(20)	
C--- user terminate routine (hztool)
cc	CALL hzraend
	
      CALL PTIME(' rapgap  ',2,0)
      CALL PTIME('        ',2,99)
      STOP
      END
      subroutine dummysub(a)
      return
      end
