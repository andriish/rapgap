************************************************************************
      SUBROUTINE FRAGPX (IBP,IDIFP)

      IMPLICIT NONE
************************************************************************
*
*     Purpose:    Take hadronic mass at p vertex and let it decay
*     --------    via an isotropic phase space decay
*
*     Input:      IBP:  Number of beam proton in /PYJETS/
*     ------      IDIFP:  Number of diffractive state in /PYJETS/
*                 +PYJETS common block
*
*
*     Called by:  PRODIFFVM
*     ----------
*
*     Author:     Benno List, 12.1.92
*     -------
*
*     Changed by: Paul Thompson
*     -----------
*
*
*     Calling:  From H1UTIL:
*     --------  RGRNG        normally distributed random number
*
*               From CERN Library:
*               V111 BINOMI  binomially distributed random numbers
*               Copied from CERN library
*
*               rgrambo        RAndom Momentum BOoster(phase space gen.)
*               RKNO:        Mean pion number given by KNO scaling
*               DECPI0       pi0 decay
*
************************************************************************
#include "rglujets.inc"
#include "rgpdiss.inc"
      INTEGER IBP,IDIFP,ISCN,NSUM,NMAX,NPI,NPM,IERROR
      INTEGER N0,NP,NM,I,J,IOK,NPOLD
      DOUBLE PRECISION DMP,DMPI,DMPI0,DMN
      REAL RNPI,RKNO,R
      Double Precision draprn
      EXTERNAL draprn
	Real Rgrng
	External Rgrng
	
* Warning
      PARAMETER (DMP = 0.938270D0,DMPI0 = 0.139568D0,DMN=1.0D0)
      PARAMETER (DMPI = 1.0D0)

      REAL*8 XM (100), PMOM (4, 100), WT
      REAL*8 PDPRS (4), PDPCMS (4)
      DOUBLE PRECISION beta(3)

*
* Determine maximal possible number of pions,
* assuming that the pions are all charged and that
* the outgoing baryon is a neutron
*

C     PRINT *, '--- FRAGPX ---'
C     PRINT *, 'NPART = ', NPART
C     PRINT *, 'PTMAXP = ', PTMAXP
     
      IF (P (IDIFP,5) .LT. DMP+DMPI0) THEN
        PRINT *, '### FRAGPX error: not enough energy! '
        PRINT *, 'IDIFP: ', IDIFP
        PRINT *, 'P (IDIFP,5): ', P(IDIFP,5)
        PRINT *, 'DMP + DMPI0:       ', DMP+DMPI0
c        CALL ERRLOG (30, 'F: FRAGPX: Not enough energy!')
        STOP
      END IF
* Determine maximal pion number

      NMAX = INT ((P(IDIFP,5) - DMN - 0.1)/DMPI)
      IF (NMAX .EQ. 0) THEN
        NMAX = 1
      ELSE IF (NMAX .GT. 99) THEN
        NMAX = 99
      END IF

c      PRINT *, 'FRAGPX: inelastic case: NMAX ', NMAX

*
* Determine multiplicity NPI (n) of pions
*

      RNPI = 1.5*RKNO (SNGL (P(IDIFP,5) - DMP))
C     PRINT *, 'mean pion number: ', RNPI

      IF (NMAX .EQ. 1) THEN
        NPI = 1
      ELSE
 20     CONTINUE
          R = RGRNG()
          R = 0.5D0*RNPI*R + RNPI + 0.5D0
	    
        IF (R .LT. 1.0 .OR. R .GE. REAL (NMAX + 1)) GOTO 20
        NPI = INT (R)
      ENDIF

c      PRINT *, 'n: ', NPI
      NSUM = NPI + 1

*
* Determine number NPM (n+-) of charged pions
* and number N0 (n0) of neutral pions so that
* about 2/3 of the pions are charged.
* Take into account that the pi+ production threshold is above
* the pi0 threshold.
* If only one pion is produced, it's 50% pi0, 50% pi+
*

      IF (NPI .GE. 2) THEN
        CALL BINOMI (NPI, 0.6666667, NPM, IERROR)
      ELSE IF (P(IDIFP,5) .LE. DMN+DMPI) THEN
        NPM = 0
      ELSE IF (draprn () .LT. 0.5) THEN
        NPM = NPI
      ELSE
        NPM = 0
      END IF
      N0 = NPI - NPM

*
* Determine the number of positive and negative pions so that
* there are no more negative than positive pions and at most
* one positive pion more than there are negative pions.
*

      NP = INT (0.5*(NPM + 1.001))
      NM = INT (0.5*(NPM + 0.999))
c      PRINT *, 'n+- = ', NPM, ' n0 = ', N0, 'n+ = ', NP, ' n- = ', NM

*
* Determine type and mass of baryon:
* if n+ = n- (pion system is neutral), outgoing baryon is a proton,
* otherwise (pion system is positive), outgoing particle is a neutron
*

      IF (NP .EQ. NM) THEN
        XM (1) = DMP
      ELSE
        XM (1) = DMN
      ENDIF

*
* Fill the mass array XM and the particle type array ITYPE
*

      DO 30, I = 2, NPM + 1
        XM (I) = DMPI
 30   CONTINUE
      DO 40, I = NPM + 2, NSUM
        XM (I) = DMPI0
 40   CONTINUE

*
* Perform isotropic phase space decay via rgrambo
*
       do 15 i=1,3
           beta(i)=P(idifp,i)/P(idifp,4)       
 15     continue
c      write(6,*) ' fragpx: before boost' 
c	call pylist(1)

       CALL PYROBO(idifp,idifp,0.d0,0.d0,-beta(1),-beta(2) 
     &                                ,-beta(3))

c      write(6,*) ' fragpx: after boost' 
c	call pylist(1)
      CALL rgrambo (NSUM, P(IDIFP,5), XM, PMOM, WT, 1, IOK)

      IF (IOK .GT. 0. and.IOK.ne.5 ) THEN
        PRINT *,'FRAGPX: rgrambo failed! MXP =', P (IDIFP,5),
     +           ', NPI =', NPI, ', IOK = ', IOK,', NSUM =', NSUM 
c        CALL ERRLOG (31, 'W: FRAGPX: rgrambo failed!')
      END IF

*
* Update /PYJETS/
*

      ISCN = N + 1

      K (IDIFP,4) = N + 1
      K (IDIFP,5) = N + NSUM
      K (IDIFP,1) = 11

      IF (NP .EQ. NM) THEN
        K (ISCN,2) = SIGN (2212, IBP)
      ELSE
        K (ISCN,2) = SIGN (2112, IBP)
      ENDIF
      DO 70, I = N + 2, N + NP + 1
        K (I,2) = SIGN (211, IBP)
 70   CONTINUE
      DO 80, I = N + NP + 2, N + NPM + 1
        K (I,2) = -SIGN (211, IBP)
 80   CONTINUE
      DO 90, I = N + NPM + 2, N + NSUM
        K (I,2) = 111
 90   CONTINUE

      DO 120, I = 1, NSUM
        P(N+i,5)=0
        DO 100, J = 1, 4
          P(N + I,J) = PMOM (J, I)
          P(N+i,5)=XM(i)
 100    CONTINUE
        K (N + I,1) = 1
        K (N + I,3) = IDIFP
        K (N + I,4) = 0
        K (N + I,5) = 0
 120  CONTINUE
       CALL PYROBO(idifp,idifp,0.d0,0.d0,beta(1),beta(2) 
     &                                ,beta(3))
c      write(6,*) ' fragpx: after boost back',nsum,n
c	call pylist(1)

*
* Perform pi0 decay
*

      NPOLD = N
      N = N + NSUM
       CALL PYROBO(npold+1,n,0.d0,0.d0,beta(1),beta(2) 
     &                                ,beta(3))

c      write(6,*) ' fragpx: after boost final boost' ,n,npold
c	call pylist(1)
      DO 130, I = NPM + 2, NSUM
       if(ijetset.eq.1)  CALL DECPI0 (NPOLD + I)
 130  CONTINUE

      RETURN
      END
