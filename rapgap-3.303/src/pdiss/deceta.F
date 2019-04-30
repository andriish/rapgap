************************************************************************
      SUBROUTINE DECETA (IETA)

      IMPLICIT NONE
************************************************************************
*
*     Purpose:    Let a eta0 decay
*     --------    eta0 -> 2 gamma (39.1%), 3 pi0 (32.2%),
*                           pi+ pi- pi0 (23.8%), pi+ pi- gamma (4.9%)
*
*     Input:      IETA:    Number of eta0 in /CPROD/
*     ------      from common block /CPROD/:
*                 NPART, PPCMS8, ITYPE
*
*     Output:     to common block /CPROD/:
*     -------     NPART, PPCMS8, ITYPE, ISTAT, MOHEP, IDAHEP
*
*     Called by:  DECPHI
*     ----------
*
*     Author:     Benno List, 7.1.94, 25.1.93 (DECAYVM)
*     -------
*
*     Changed by:
*     -----------
*
*     Calling:  
*               draprn         Random number generator
*
*               rgrambo        RAndom Momenta BOoster (phase space gen.)
*               DECPI0       pi0 decay
*
************************************************************************
#include "rglujets.inc"
#include "rgludat2.inc"
#include "rgpdiss.inc"
	
      Double Precision xm1,xw1
      INTEGER IETA,N1,N2,ITER,J,IOK,NSUM,I
      REAL R
      Double Precision draprn
      EXTERNAL draprn
      REAL*8 XM (100), PMOM (4, 100), WT, SUM
      REAL*8 PDPCMS (1:4), PDPRS (1:4)
      LOGICAL ACCEPT

      DOUBLE PRECISION beta(3)

      if(ijetset.eq.1) return
C     PRINT *, '--- DECETA ----'
C     CALL PR5V8 ('p_eta0', PPCMS8 (1, IETA))

      IF (K (IETA,2) .NE. 221) THEN
        PRINT *, 'DECETA error: ITYPE (', IETA, ') = ', K (IETA,2)
c        CALL ERRLOG (190, 'S: DECETA: Wrong particle!')
        RETURN
      END IF
      IF (K (IETA,1) .NE. 1) THEN
        PRINT *, 'DECETA error: ISTAT (', IETA, ') = ', K(IETA,1)
c        CALL ERRLOG (191, 'W: DECETA: eta not undecayed!')
        RETURN
      END IF

*
* Choose decay channel
*

      R = draprn ()
      IF (R .LT. 0.049) THEN
        NSUM = 3
        K (N+1,2) =  211
        K (N+2,2) = -211
        K (N+3,2) =   22
      ELSE IF (R .LT. 0.287) THEN
        NSUM = 3
        K (N+1,2) =  211
        K (N+2,2) = -211
        K (N+3,2) =  111
      ELSE IF (R .LT. 0.287) THEN
        NSUM = 3
        K (N+1,2) =  111
        K (N+2,2) =  111
        K (N+3,2) =  111
      ELSE
        NSUM = 2
        K (N+1,2) =   22
        K (N+2,2) =   22
      END IF

*
* Fill the mass array XM
*

      SUM = 0D0
      DO 10, I = 1, NSUM
	  call pdismass(K(N+I,2),XM1,XW1)
        XM (I) = XM1
        SUM = SUM + XM (I)
 10   CONTINUE

      IF (SUM .GT. P (IETA,5)) THEN
        PRINT *, 'DECETA error: mass ', P(IETA,5), 'too small!'
c        CALL ERRLOG (192, 'S: DECETA: mass of eta too small!')
        RETURN
      END IF

*
* Iterate rgrambo call until cos (theta*) fits distribution;
* maximally 100 iterations are allowed.
*
      
       do 15 i=1,3
           beta(i)=P(ieta,i)/P(ieta,4)       
 15     continue

       CALL pyrobo(ieta,ieta,0.d0,0.d0,-beta(1),-beta(2) 
     &                                ,-beta(3))

      ITER = 0
      
 20   CONTINUE
      ITER = ITER + 1

*
* Call rgrambo: distribute particles uniformely in phase space
*

      CALL rgrambo (NSUM, P (IETA,5), XM, PMOM, WT, 1, IOK)
      IF (IOK .NE. 0) THEN
c         CALL ERRLOG (193, 'S: DECETA: rgrambo failed!')
      END IF
*
* Calculate cos (theta*)
*

      IF (NSUM .EQ. 2) THEN
         ACCEPT = .TRUE.
      ELSE
         ACCEPT = .TRUE.
      END IF

      IF (.NOT. ACCEPT .AND. ITER .LE. 100) GO TO 20



      IF (ITER .GT. 100) THEN
         PRINT *, '*** DECETA warning: more than 100 rgrambo iterations'
c         CALL ERRLOG (194, 'S: DECETA: mass of eta too small!')
      END IF

      PRINT *, 'deceta: Iterations: ', ITER

*
* Update /DUJETS/
*

       do 30 i=1,nsum            
          P(N+i,5)=0
          do 40 j=1,4      
             P(N+i,j)=PMOM(j,i)
c             P(N+i,5)=P(N+i,5)+PMOM(j,i)
             P(N+i,5)=xm(i)
 40       continue
c          P(N+i,4)=PMOM(4,i)
c          P(N+i,5)=P(N+i,4)-P(N+i,5)

 30    continue 


       CALL pyrobo(ieta,ieta,0.d0,0.d0,beta(1),beta(2) 
     &                                ,beta(3))

       CALL pyrobo(N+1,N+NSUM,0.d0,0.d0,beta(1),beta(2) 
     &                                ,beta(3))



      write(*,*) 'rgrambo eta'
      call pylist(1)


      K (IETA,1) = 11
      K (IETA,4) = N + 1
      K (IETA,5) = N + NSUM

      DO 50, I = 1, NSUM
*        P (N + I,5) = XM (I)
        K(N + I,1) = 1
        K(N + I,3) = IETA
        K(N + I,4) = 0
        K(N + I,5) = 0
 50   CONTINUE
      N = N + NSUM

*
* Perform pi0 decay
*
      N1 = N - NSUM + 1
      N2 = N
      DO 60, I = N1, N2
         IF (K(I,2) .EQ. 111) THEN
            CALL DECPI0 (I)
         END IF
 60   CONTINUE

      K (IETA,1) = 11 

      RETURN
      END
