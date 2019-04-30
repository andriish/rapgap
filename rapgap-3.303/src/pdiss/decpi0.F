************************************************************************
      SUBROUTINE DECPI0 (IPI0)

      IMPLICIT NONE
************************************************************************
*
*     Purpose:   Let a pi0 decay isotropically into photons
*     --------
*
*     Input:     IPI0:   Number of pi0 in /PYJETS/
*                
*     ------     
*
*     Output:    to common block /PYJETS/
*     -------    
*
*     Author:    Benno List  3.2.94
*     -------
*
*     Called by: 
*     ----------
*
*
************************************************************************
#include "rglujets.inc"
#include "rgpdiss.inc"
************************************************************************
      INTEGER IPI0
      REAL*8 PDPRS (4), CTHETA, PHI, E
      INTEGER ISIG,I,J,IOK,NOLD
      DOUBLE PRECISION beta(3)

       REAL*8 PMOM (4, 100), WT, SUM,xm(100)
      if(ijetset.eq.0) return

      IF (K (IPI0,2) .NE. 111) THEN
        PRINT *, 'DECPI0 error: ITYPE (', IPI0, ') = ', K (IPI0,2)
        RETURN
      END IF
      IF (K (IPI0,1) .NE. 1) THEN
        PRINT *, 'DECPI0 error: ISTAT (', IPI0, ') = ', K (IPI0,1)
        RETURN
      END IF
      K (IPI0,1) = 11
      K (IPI0,4) = N+1
      K (IPI0,5) = N+2

*******************************

* Call rgrambo: distribute particles uniformely in phase space
*

       do 15 i=1,3
           beta(i)=P(ipi0,i)/P(ipi0,4)       
 15     continue

       CALL PYROBO(ipi0,ipi0,0.d0,0.d0,-beta(1),-beta(2),-beta(3))


c      write(6,*) ' decpi0 ',ipi0,P(Ipi0,5)
c	call pylist(1)
	xm(1) = 0d0
	xm(2) = 0d0
      CALL rgrambo (2, P(Ipi0,5), XM, PMOM, WT, 1, IOK)

      IF (IOK .NE. 0) THEN
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
c       write(6,*) ' before boost '
c       write(*,*) 'decpi0',P(N+1,1),P(N+1,2),P(N+1,3),P(N+1,4),P(N+1,5) 
c       write(*,*) 'decpi0',P(N+2,1),P(N+2,2),P(N+2,3),P(N+2,4),P(N+2,5)

       CALL PYROBO(ipi0,ipi0,0.d0,0.d0,beta(1),beta(2),beta(3))
     
c       write(6,*) ' after boost '
c       write(*,*) 'decpi0',P(N+1,1),P(N+1,2),P(N+1,3),P(N+1,4),P(N+1,5) 
c       write(*,*) 'decpi0',P(N+2,1),P(N+2,2),P(N+2,3),P(N+2,4),P(N+2,5)
      NOLD = N 
      N = N+1
*      P(N,5) = 0.135
      K(N,1) = 1
      K(N,2) = 22
      K(N,3) = IPI0
      K(N,4) = 0
      K(N,5) = 0

      N = N+1
*     P (N,5) = 0.135
      K(N,1) = 1
      K(N,2) = 22
      K(N,3) = IPI0
      K(N,4) = 0
      K(N,5) = 0
      CALL PYROBO(NOLD+1,N,0.d0,0.d0,beta(1),beta(2),beta(3))
c      write(*,*) 'rgrambo pi0'
c      call pylist(1)

      RETURN
      END


