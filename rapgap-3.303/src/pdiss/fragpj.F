************************************************************************
      SUBROUTINE FRAGPJ (IBP, IDIFP)

      IMPLICIT NONE
************************************************************************
*
*     Purpose:    Take hadronic mass at p vertex and let it decay
*     --------    
*
*     Input:      IBP:    Number of beam proton in PYJETS
*     ------      IDIFP:  Number of diffractive state in PYJETS
*
*     Called by:  PRODIFFVM
*     ----------
*
*     Author:     Benno List, 17.1.94
*     -------
*
*     Changed by: Paul Thompson
*     -----------
*
*     Calling:    SPLITP
*     --------    
*
************************************************************************
#include "rglujets.inc"
#include "rgpdiss.inc"
********************************************************
      DOUBLE PRECISION DMP,DMPI0
      PARAMETER (DMP = 0.938270D0,DMPI0 = 0.139568D0)
      DOUBLE PRECISION PQQ (5), PQ (5), PPOM (5)
      INTEGER IBP, IDIFP
********************************************************
      DOUBLE PRECISION PPCMS8(5)
      Double Precision draprn
      EXTERNAL draprn
      INTEGER IPOM,IBEAMP,IQ,IQQ,I,IERR,ICHECK
      REAL R
      DATA ICHECK /0/

      Integer npdc1,npdc2
      common/pdcyn/npdc1,npdc2

      DOUBLE PRECISION la_pg(5),la_prem(5)
	DOUBLE PRECISION ptes(5)




*****
      INTEGER IARR12 (2) /1, 2/
      
*****


c     Integer I_glu
c	Data I_glu/0/  !  old diffvm treatment 
c	Data I_glu/1/  !  new treatment, adding a gluon, done by Leszek Adamczyk



********************************************************

C     PRINT *, '--- FRAGPJ ---'
C     PRINT *, 'IBP   = ', IBP
C     PRINT *, 'IDIFP = ', IDIFP

c      write(6,*) ' fragjp I_glu :',I_glu

      IF ( P(IDIFP,5) .LT. DMP+DMPI0) THEN
         PRINT *, '### FRAGPJ error: not enough energy! '
         PRINT *, 'IBP:   ', IBP
         PRINT *, 'IDIFP: ', IDIFP
         PRINT *, 'P(IDIFP,5): ', P(IDIFP,5)
         PRINT *, 'DMP + DMPI0:       ', DMP+DMPI0
c        CALL ERRLOG (40, 'F: FRAGPJ: Not enough energy!')
      END IF
*
* Split proton into q - qq system: choose q and qq codes
* The pomeron couples with the same strength to u and d quarks.
* ud_1 and uu_1 are spin triplets, ud_0 is a spin singulett.
* Therefore q-qq is 1/7 u-ud_0, 3/7 u-ud_1, and 3/7 d-uu_1.
*
      IBEAMP=K(IBP,2)
      IF (ICHECK.EQ.1) write(6,*) 'SPLITP: ibeamp',ibeamp
      R = 7.0*draprn()
      IF (R .LE. 1.0) THEN
         IQ  = SIGN (2, IBEAMP)       !u
         IQQ = SIGN (2101, IBEAMP)    !ud spin 0
      ELSE IF (R .LE. 4.0) THEN
         IQ  = SIGN (2, IBEAMP)       !u
         IQQ = SIGN (2103, IBEAMP)    !ud spin 1
      ELSE
         IQ  = SIGN (1, IBEAMP)       !d
         IQQ = SIGN (2203, IBEAMP)    !uu spin 1
      END IF
      
      IF (ICHECK.EQ.1) THEN
      write (6,*) 'FRAGPJ IQ:  ', IQ,  ', mass: ', PYMASS (IQ)
      write (6,*) 'FRAGPJ IQQ: ', IQQ, ', mass: ', PYMASS (IQQ)
      ENDIF

      DO I=3,N
         IF(IABS(K(I,2)).EQ.990) IPOM = I
      enddo
* Fill 5 vector of the beam proton and pomeron 
      do i=1,5
         PPCMS8(I)= P(IBP,I)
chju         PPOM (I) = P(IPOM,I)
* To avoid numerical problems, calculate pomeron momentum from
* proton side
	   PPom(I) = P(IDIFP,I) - P(IBP,I)
         IF (ICHECK.EQ.1) THEN
         write(6,*)'I,PPOM,P(IPOM),PROTON',I,PPOM(I),P(IPOM,I),P(IBP,I)
         ENDIF
      enddo
      PPOM (5) = P(IPOM,5)

      IF (ICHECK.EQ.1) THEN
         write(6,*) 'FRAGPJ: PPOM(5)',PPOM(5)
         write(6,*) 'FRAGPJ: CALL SPLITP'
      ENDIF

      If(I_glu.eq.0) then 
          CALL SPLITP (PPCMS8, PPOM, PQ, PQQ, IERR)
      else
cnew
          CALL SPLITP (PPCMS8, PPOM, la_pg, la_prem, IERR)
          do i=1,5
           pq(i) = 1./3. * la_prem(i)
           pqq(i) = 2./3. * la_prem(i)
	     ptes(i)=pq(i)+pqq(i)+la_pg(i)
          enddo
c	    write(6,*) 'ptest ',ptes
	endif
cnew
      IF (ICHECK.EQ.1) THEN
         write (6,*) 'N before',N
         write (6,*) 'SPLITP',IERR
      ENDIF

* Proton  decay

 
c      write(*,*) 'pylist 1',Ierr
c      call pylist(1)

c    write(6,*) 'PRODIFFVM: idifp, ipom, N', idifp,ipom,n

      K(IDIFP,1)=11
      N=N+1
      npdc1 = N
      DO I=1,5
         P(N,I)=PQQ (I)
      ENDDO
      K(N,1)=12
cy      K(N,1)=11  ! diquark 
      K(N,2)=IQQ
      K(N,3)=IDIFP

c      IF (ICHECK.EQ.1) write(6,*)'in between N',N
      If(I_glu.eq.1) then
cnew
         N=N+1
         DO I=1,5
            P(N,I)=la_Pg(I)
         ENDDO
         K(N,1)=12
         K(N,2)=21
         K(N,3)=IDIFP
cnew	
      Endif
      N=N+1
      npdc2 = N
      DO I=1,5
         P(N,I)=PQ(I)
      ENDDO
      K(N,1)=11
cy      K(N,1)=11 ! quark
      K(N,2)=IQ
      K(N,3)=IDIFP



      IF (ICHECK.EQ.1) THEN
         write(6,*)'after adding 2'
         write(6,*)'N',N
         write(6,*) 'Fragpj: end pylist'
         CALL pylist (1)
      ENDIF

c      write(*,*) 'FRAGPJ'
c      call pylist(1)
   
      RETURN
      END





