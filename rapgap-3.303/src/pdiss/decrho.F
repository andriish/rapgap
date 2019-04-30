************************************************************************
      SUBROUTINE DECRHO (IRHO, HELI)

      IMPLICIT NONE
************************************************************************
*
*     Purpose:    Let a rho decay
*     --------    rho+ -> pi+ pi0
*                 rho- -> pi- pi0
*                 rho0 -> pi+ pi- (99.0%), pi+ pi- gamma (1.0%)
*
*     Input:      IRHO:    Number of rho in /CPROD/
*     ------      HELI:    Helicity of the rho meson
*
*     Called by:  
*     ----------
*
*     Author:     Benno List, 7.1.94, 25.1.93 (DECAYVM)
*     -------
*
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
      INTEGER IRHO, HELI,J,ITER,NSUM,IPI0,I,IOK,NOLD

      DOUBLE PRECISION boost(3)
      Double Precision draprn
      EXTERNAL draprn

      REAL CTHST,R,ALPHA,BETA,PVM2
      REAL*8 XM (100), PMOM (4, 100), WT, SUM
      REAL*8 PDPCMS (1:4), PDPRS (1:4)
      LOGICAL ACCEPT

C     PRINT *, '--- DECRHO ----'
C     CALL PR5V8 ('p_rho ', PPCMS8 (1, IRHO))
      if(ijetset.eq.0) return

      IF (ABS (K(IRHO,2)) .NE. 213 .AND. K (IRHO,2) .NE. 113) THEN
        PRINT *, 'DECRHO error: ITYPE (', IRHO, ') = ', K (IRHO,2)
c        CALL ERRLOG (160, 'S: DECRHO: Wrong particle!')
        RETURN
      END IF
      IF (K(IRHO,1) .NE. 1) THEN
        PRINT *, 'DECRHO error: ISTAT (', IRHO, ') = ', K(IRHO,1)
c        CALL ERRLOG (161, 'W: DECRHO: rho not undecayed!')
        RETURN
      END IF

*
* Determine parameters for decay angle distribution
*

      IF (HELI .EQ. 0) THEN
        ALPHA = 0.0
        BETA  = 1.0
      ELSE
        ALPHA = 1.0
        BETA  = 0.0
      END IF

*
* Choose decay channel
*

      R = draprn ()
      IPI0 = 0
      IF (K(IRHO,2) .EQ. 113) THEN
        IF (R .LT. 0.010) THEN
          NSUM = 3
          K (N+1,2) =  211
          K (N+2,2) = -211
          K (N+3,2) =   22
        ELSE
          NSUM = 2
          K (N+1,2) =  211
          K (N+2,2) = -211
        END IF
      ELSE
        NSUM = 2
        K (N+1,2) =  SIGN (211, K (IRHO,2))
        K (N+2,2) =  111
        IPI0 = N + 2
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

      IF (SUM .GT. P (IRHO,5)) THEN
        PRINT *, 'DECRHO error: mass ', P(IRHO,5), 'too small!'
c        CALL ERRLOG (162, 'S: DECRHO: mass of rho too small!')
        RETURN
      END IF

*
* Iterate rgrambo call until cos (theta*) fits distribution;
* maximally 100 iterations are allowed.
*

      ITER = 0


        do 15 i=1,3
           boost(i)=P(irho,i)/P(irho,4)       
 15     continue

       CALL pyrobo(irho,irho,0.d0,0.d0,-boost(1),-boost(2),-boost(3))
      PVM2 = P (IRHO,1)**2 + P (IRHO,2)**2 +
     +     P (IRHO,3)**2

 20   CONTINUE
      ITER = ITER + 1

*
* Call rgrambo: distribute particles uniformely in phase space
*


       CALL rgrambo (NSUM, P(IRHO,5), XM, PMOM, WT, 1, IOK)
       IF (IOK .NE. 0) THEN
         write(6,*) ' DECNST: rgrambo failed!'
       END IF


*
* Calculate cos (theta*)
*

        IF (NSUM .EQ. 2) THEN
          CTHST = SNGL (P (IRHO,1)*PMOM (1, 1) +
     +                  P (IRHO,2)*PMOM (2, 1) +
     +                  P (IRHO,3)*PMOM (3, 1)) /
     +          SQRT (SNGL ((PMOM (1, 1)**2+PMOM (2, 1)**2
     +          +PMOM (3, 1)**2)*PVM2))

          IF (CTHST .LT. -1.0) THEN
            PRINT *, '*** DECRHO warning: CTHST < -1: ', CTHST
C            CALL ERRLOG (164, 'W: DECRHO: CTHST < -1!')
          ELSE IF (CTHST .GT. 1.0) THEN
            PRINT *, '*** DECRHO warning: CTHST > 1: ', CTHST
C            CALL ERRLOG (165, 'W: DECRHO: CTHST > -1!')
          END IF
          ACCEPT = (draprn () .LT. ALPHA + (BETA-ALPHA)*CTHST**2)
c           write(6,*) ALPHA + (BETA-ALPHA)*CTHST**2,CTHST
        ELSE
          ACCEPT = .TRUE.
        END IF

      IF (.NOT. ACCEPT .AND. ITER .LE. 100) GO TO 20

      IF (ITER .GT. 100) THEN
        PRINT *, '*** DECRHO warning: more than 100 rgrambo iterations'
c        CALL ERRLOG (166, 'W: DECRHO: More than 200 rgrambo iterations!')
      END IF

C     PRINT *, 'Iterations: ', ITER

*
* Update /PYJETS/
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
 
c       write(*,*) ' decrho: p(instar) ', P(irho,1),P(irho,2),
c     &              P(irho,3),P(irho,4),P(irho,5)
c       do i=1,nsum
c       write(*,*) ' decrho: p(n+i) ',i,P(N+i,1),P(N+i,2),P(N+i,3),
c     & P(N+i,4),P(N+i,5) 
c       enddo

       CALL pyrobo(IRHO,IRHO,0.d0,0.d0,boost(1),boost(2),boost(3))

       NOLD = N
c      write(*,*) 'rgrambo rho'
c      call pylist(1)



      K (IRHO,1) = 2
      K (IRHO,4) = N + 1
      K (IRHO,5) = N + NSUM

      DO 50, I = 1, NSUM
*        P(N+I,5) = XM (I)
        K (N + I,1) = 1
        K (N + I,3) = IRHO
        K (N + I,4) = 0
        K (N + I,5) = 0
 50   CONTINUE

      N = N + NSUM
      CALL pyrobo(NOLD+1,N,0.d0,0.d0,boost(1),boost(2),boost(3))

      K (IRHO,1) = 11
c      write(*,*) 'decrho  rho'
c      call pylist(1)

*
* Perform pi0 decay
*

      IF (IPI0 .GT. 0) THEN
        CALL DECPI0 (IPI0)
      END IF

             
      RETURN
      END
