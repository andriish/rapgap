      SUBROUTINE PRODIFFVM

      IMPLICIT NONE
*
*  PROTON DISSOCIATION A LA DIFFVM
*
#include "rglujets.inc"
#include "rgpdiss.inc"
#include "rgpart.inc"

      DOUBLE PRECISION DMXP       
      INTEGER IBP,IDIFP,J,ICHECK,IPROT,I
      INTEGER IPOM,KPART
      REAL RFRAC
      Double Precision draprn
      EXTERNAL draprn
	
      DATA ICHECK /0/


      Integer npdc1,npdc2
      common/pdcyn/npdc1,npdc2  

      DOUBLE PRECISION CMGP(4)

      integer pdvm
      common/prodvm/pdvm

cy in ep CMS FRAME

* Type of dissociation
c      IFRAGP=2

* beam proton
      DO 5 I=1,2
         IF(IABS(K(I,2)).EQ.2212) IBP = I
 5    CONTINUE
 
      IF (ICHECK.EQ.1) THEN
         write(6,*)'PRODIFFVM: beam proton ibp',IBP
      ENDIF
 
      DO 10 I=3,N
         IF(IABS(K(I,2)).EQ.2212) IPROT = I
         IF(IABS(K(I,2)).EQ.990) IPOM = I
 10   CONTINUE
      IF(IPROT.LT.0) RETURN
      IDIFP=IPROT
      KPART = K(IPROT,2)  

      DMXP=P(IDIFP,5)


c      write(6,*) 'prodiffvm MY = ', DMXP, sqrt(MY2),ifragp

* Determine the state
      IF (IFRAGP .EQ. 1 .OR. IFRAGP .EQ. 2) THEN
         IF (DMXP .LT. 1.48) THEN
            K(IDIFP,2) = SIGN (12212,KPART)
            K(IDIFP,1)=11
            npdc1 = idifp
            npdc2 = 0
         ELSE IF (DMXP .LT. 1.60) THEN
            K(IDIFP,2) = SIGN (2124,KPART)
            K(IDIFP,1)=11
            npdc1 = idifp
            npdc2 = 0
         ELSE IF (DMXP .LT. 1.90) THEN
            RFRAC = draprn()
            IF (RFRAC .LT. 0.5) THEN
               K(IDIFP,2)   = SIGN (12216, KPART)
               K(IDIFP,1)=11
               npdc1 = idifp
               npdc2 = 0
            ELSE IF (RFRAC .LT. 0.83) THEN
               K(IDIFP,2) = SIGN (22124, KPART)
               K(IDIFP,1)=11
                npdc1 = idifp
                npdc2 = 0
            ELSE
               K(IDIFP,2) = SIGN (42212, KPART)
               K(IDIFP,1)=11
               npdc1 = idifp
               npdc2 = 0
            END IF
         ELSE
* make proton a diffracted p_diffr+ state 
c            K(IDIFP,2) = SIGN (2210, KPART)
            K(IDIFP,2) = SIGN (9902210, KPART)
            pdvm=1
         END IF
      ELSE
         K(IDIFP,2) = SIGN (ABS (IFRAGP), KPART)
         write(6,*) 'PRODIFFVM:   incorrect value of IFRAGP'
      END IF

      DO 20 I=1,4
         CMGP(I) =P(2,I)+P(3,I)
   20 CONTINUE   

* THE FOLLOWING IS FROM GENEVT IN DIFFVM
*
* If necessary, perform fragmentation of diffractive proton system
*
      IF (ABS (K (IDIFP,2)) .EQ. 2212) THEN
         write(6,*) 'PRODIFFVM:Should not be true:Elastic proton!'
      ELSE IF (ABS (K (IDIFP,2)) .EQ. 9902210) THEN
         IF (IFRAGP .EQ. 1) THEN   
c	      write(6,*) ' prodiffvm: no frag. for FRAGPJ '         
            CALL FRAGPJ (IBP,IDIFP)
         ELSE
c	      write(6,*) ' prodiffvm: no frag. for FRAGPX '         
            CALL FRAGPX (IBP,IDIFP)
         END IF
      endif

* do the decay of the resonances later on in pdisdcvm.f


      RETURN
      END




