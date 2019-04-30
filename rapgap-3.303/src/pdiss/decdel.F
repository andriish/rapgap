************************************************************************
      SUBROUTINE DECDEL (IDELTA)

      IMPLICIT NONE
************************************************************************
*
*     Purpose:    Let a Delta(1232) decay
*     --------    Delta++ -> p pi+ (100%)
*                 Delta+  -> p pi0 (67%), n pi+ (33%)
*                 Delta0  -> n pi0 (67%), p pi- (33%)
*                 Delta-  -> n pi- (100%)
*
*     Input:      IDELTA:  Number of Delta in /PYJETS/
*     ------      
*
*     Output:    
*     -------     
*
*     Called by:  DECNST
*     ----------
*
*     Author:     Benno List, 1.9.94
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
**************************************************
      INTEGER IDELTA,ISIG,I,J,IOK,NOLD
      Double Precision xm1,xw1
	Double Precision draprn
      EXTERNAL draprn

      REAL*8 XM (100), PMOM (4, 100), WT, SUM
      REAL*8 PDPCMS (1:4), PDPRS (1:4)

      DOUBLE PRECISION beta(3)

C     PRINT *, '--- DECDEL ----'
C     CALL PR5V8 ('p_Del ', PPCMS8 (1, IDELTA))
C     PRINT *, ITYPE (IDELTA), ' decays to'

      IF (K (IDELTA,1) .NE. 1) THEN
        PRINT *, 'DECDEL error: ISTAT (', IDELTA, ') = ', K(IDELTA,1)
c        CALL ERRLOG (210, 'W: DECDEL: Particle not undecayed!')
        
        RETURN
      END IF

*
* Fill ITYPE
*
      ISIG = SIGN (1, K(IDELTA,2))
C     PRINT *, 'ISIG = ', ISIG
      IF (ABS (K(IDELTA,2)) .EQ. 2224) THEN
         K (N+1,2) = 2212*ISIG
         K (N+2,2) =  211
      ELSE IF (ABS (K (IDELTA,2)) .EQ. 2214) THEN
         IF (draprn () .LE. 0.333) THEN
            K (N+1,2) = 2112*ISIG
            K (N+2,2) =  211*ISIG
         ELSE
            K (N+1,2) = 2212*ISIG
            K (N+2,2) =  111
         END IF
      ELSE IF (ABS (K (IDELTA,2)) .EQ. 2114) THEN
         IF (draprn () .LE. 0.333) THEN
            K (N+1,2) = 2212*ISIG
            K (N+2,2) = -211*ISIG
         ELSE
            K (N+1,2) = 2112*ISIG
            K (N+2,2) =  111
         END IF
      ELSE IF (ABS (K (IDELTA,2)) .EQ. 1114) THEN
         K (N+1,2) = 2112*ISIG
         K (N+2,2) = -211*ISIG
      ELSE
         
         PRINT *, '*** DECDEL Error: ', K (IDELTA,2), ' is not a Delta'
c         CALL ERRLOG (211, 'S: DECDEL: Wrong particle!')
         RETURN
      END IF

*
* Fill the mass array XM
*

      SUM = 0D0
      DO 10, I = 1, 2
	   call pdismass(K(N+I,2),XM1,XW1)
c        XM (I) = DBLE (PMASS (K (N+I,2)))
        XM (I) = XM1
        SUM = SUM + XM (I)
 10   CONTINUE


      IF (SUM .GT. P(IDELTA,5)) THEN
        PRINT *, 'DECDEL error: mass ', P (IDELTA,5), 'too small!'
c        CALL ERRLOG (212, 'S: DECDEL: Mass of Delta too small!')
        RETURN
      END IF

*
* Call rgrambo: distribute particles uniformely in phase space
c      write(*,*) ' decdel: before boosts   '
c      call pylist(1)
*

       do 15 i=1,3
           beta(i)=P(idelta,i)/P(idelta,4)       
 15     continue

       CALL PYROBO(idelta,idelta,0.d0,0.d0,-beta(1),-beta(2),-beta(3))


c      write(*,*) ' decdel: after boost delta to cms    '
c      call pylist(1)

      CALL rgrambo (2, P(IDELTA,5), XM, PMOM, WT, 1, IOK)

      IF (IOK .NE. 0) THEN
c        CALL ERRLOG (213, 'S: DECDEL: rgrambo failed!')
        write(6,*) 'rgrambo FAILED'
      END IF



       do 20 i=1,2            
          P(N+i,5)=0
          do 40 j=1,4      
             P(N+i,j)=PMOM(j,i)
c             P(N+i,5)=P(N+i,5)+PMOM(j,i)
             P(N+i,5)=xm(i)
 40       continue
c          P(N+i,4)=PMOM(4,i)
c          P(N+i,5)=P(N+i,4)-P(N+i,5)

 20    continue   
 

c       write(*,*) ' decdel: p(instar) ', P(instar,1),P(instar,2),
c     &              P(instar,3),P(instar,4),P(instar,5)
c       do i=1,nsum
c       write(*,*) ' decdel: p(n+i) ',i,P(N+i,1),P(N+i,2),P(N+i,3),
c     & P(N+i,4),P(N+i,5) 
c       enddo

     
       CALL PYROBO(idelta,idelta,0.d0,0.d0,beta(1),beta(2),beta(3))



*
* Update /PYJETS/
*
      NOLD=N
      K(IDELTA,1) = 11
      K(IDELTA,4) = N + 1
      K(IDELTA,5) = N + 2


      DO 50, I = 1, 2
        P (N+I,5) = XM (I)
        K (N+I,1) = 1
        K (N+I,3) = IDELTA
        K (N+I,4) = 0
c        K (N+I,5) = 0
 50   CONTINUE
      N = N + 2
c      write(*,*) 'decdel after  delta boosts  ',N
c      call pylist(1)
       CALL PYROBO(NOLD+1,N,0.d0,0.d0,beta(1),beta(2),beta(3))

c      write(*,*) 'decdel after  all boosts  '
c      call pylist(1)

*
* Perform pi0 decay
*

      IF (K (N,2) .EQ. 111) THEN
        if(ijetset.eq.1) CALL DECPI0 (N)
      END IF


      
      RETURN
      END

