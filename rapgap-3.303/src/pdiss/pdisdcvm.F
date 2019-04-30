      SUBROUTINE PDISDCVM

      IMPLICIT NONE
*
*  PROTON DISSOCIATION A LA DIFFVM
*
#include "rglujets.inc"
#include "rgherup.inc"

      DOUBLE PRECISION DMXP       
      INTEGER IFRAGP,IBP,IDIFP,J,ICHECK,IPROT,I
      INTEGER IPOM,KPART
      REAL RFRAC
      Double Precision draprn
      EXTERNAL draprn
      DATA ICHECK /1/

      DOUBLE PRECISION SSS,CM,DBCMS
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)


      Integer npdc1,npdc2
      common/pdcyn/npdc1,npdc2  

      DOUBLE PRECISION CMGP(4)
      
      DOUBLE PRECISION beta(3)

      double precision my2
      integer pdvm
      common/mykine/ MY2
      common/prodvm/pdvm

cy in ep CMS FRAME

* Type of dissociation
      IFRAGP=1
      pdvm=0
* beam proton

      if(npdc1.eq.0.and.npdc2.eq.0) write(6,*) ' pdisdcvm ',npdc1,npdc2
      if(npdc2.ne.0) then
         K(npdc1,1) = 2
	   if(npdc1+1.ne.npdc2) then
	     K(npdc1+1,1) = 2
	   endif
         K(npdc2,1) = 1

         if(ILHA.eq.0) call pyexec

      else

* decay the resonances
c        K(npdc1,1) = 1

c        write(*,*) 'pylist in pdiscvm '
c        call pylist(1)
 
        IF (ABS (K (npdc1,2)) .EQ. 12212 .OR.
     +        ABS (K (npdc1,2)) .EQ.  2124 .OR.
     +        ABS (K (npdc1,2)) .EQ. 22212 .OR.
     +        ABS (K (npdc1,2)) .EQ. 32212 .OR.
     +        ABS (K (npdc1,2)) .EQ.  2216 .OR.
     +        ABS (K (npdc1,2)) .EQ. 12216 .OR.
     +        ABS (K (npdc1,2)) .EQ. 22124 .OR.
     +        ABS (K (npdc1,2)) .EQ. 42212 .OR.
     +        ABS (K (npdc1,2)) .EQ. 32124) THEN
c           write(6,*)  'K(IDIFP,2)=', K(IDIFP,2)
c                 write(6,*) 'pdisdcvm: no decnst '
              K(npdc1,1) = 1
              CALL DECNST (npdc1)
         ELSE IF (ABS (K (IDIFP,2)) .EQ. 2214) THEN
              K(npdc1,1) = 1
c                 write(6,*) 'pdisdcvm: no decdel '
              CALL DECDEL (npdc1)
         ELSE
         END IF
 
      endif

      RETURN
      END
