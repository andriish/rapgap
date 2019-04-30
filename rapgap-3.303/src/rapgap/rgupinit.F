C...RGUPINIT
 
      SUBROUTINE RGUPINIT
      Implicit None
#include "rglujets.inc"
#include "rgeffic.inc"
#include "rgluco.inc"
#include "rgpara.inc"
#include "rgherup.inc"
C...RAPGAP process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/
	
	Integer MSTP,MSTI
	Double Precision PARP,PARI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      save /pypars/
      Integer MINT
      Double PRecision VINT
      COMMON/PYINT1/MINT(400),VINT(400)
      LOGICAL ex
      write(6,*) ' RAPGAP rgupinit'	
C Decide on weighting strategy, events have weight = 1 
	IDWTUP = 3
C Nr of external processes
      NPRUP = 1
c set incoming beams
      IDBMUP(1) = K(2,2)
      IDBMUP(2) = K(1,2)
	
	LPRUP(1)=IPRO
	EBMUP(1) = P(2,4)+1D-7*P(2,4)
	EBMUP(2) = P(1,4)+1D-7*P(1,4)
      PDFGUP(1) = -1
      PDFGUP(2) = -1
	PDFSUP(1) =-1.
	PDFSUP(2) =-1.
	XSECUP(1) = Avgi*1000.D0
	XERRUP(1) = SD*1000.D0
	XMAXUP(1) = 1.
C write xsection information
      IF(ILHA.ge.10) then
c call directly PYTHIA to do final state radiation and fragmentation	
	   MSTP(61)=0 ! initial state PS
c         mstp(71)=1
c	   if(ifps.le.1) mstp(71)=0
c      MSTP(71)=0 ! final state PS
         MSTP(81)=0 ! multiple interactions
c      MSTP(111)=0 ! fragmentation 
	   IF(NFRAG.eq.0) MSTP(111)=0
C...Initialize.
         CALL PYINIT('USER',' ',' ',0D0)
         write(6,*) ' rgupinit ',VINT(1),VINT(2)
c         pause
        IDWTUP   = 1
        write(6,*) ' end of pythia ini ',mstp(81)
	else
         inquire(FILE='RAPGAP.gen',EXIST=ex)
         if(ex) then
	   Write(6,*) ' File RAPGAP.gen exists, will be overwritten '
         OPEN(21,FILE='RAPGAP.gen',FORM='unformatted',STATUS=
     +         'OLD',ERR=100)
	   else
	   Write(6,*) ' Open new file RAPGAP.gen '
         OPEN(21,FILE='RAPGAP.gen',FORM='unformatted',STATUS=
     +         'NEW',ERR=100)
         endif
         WRITE(21,ERR=100) IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),
     &   PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),
     &   IDWTUP,NPRUP,XSECUP(1),XERRUP(1),XMAXUP(1),
     &   LPRUP(1)
      endif 	
      RETURN
100	Write(6,*) ' Error! File open failed for RAPGAP.gen  STOP '
      STOP
      END 
