************************************************************************
      SUBROUTINE DECNST (INSTAR)

      IMPLICIT NONE
************************************************************************
*
*     Purpose:    Let a N* decay
*     --------    Possible decay modes (also charge conjugates)
*                 N*+ -> p pi0, n pi+
*                 N*+ -> Delta++ pi-, Delta+ pi0, Delta0 pi+
*                 N*+ -> p (pi+ pi-)_I=0, p (pi0 pi0)_I=0
*                 N*+ -> p rho0, n rho+
*                 N*+ -> p eta
*                 N*+ -> Lambda K+
*                 N*0 -> p pi-, n pi0
*                 N*0 -> Delta+ pi-, Delta0 pi0, Delta- pi+
*                 N*0 -> n (pi+ pi-)_I=0, n (pi0 pi0)_I=0
*                 N*0 -> p rho-, n rho0
*                 N*0 -> n eta
*                 N*0 -> Lambda K0S, Lambda K0L
*                 Branching ratios between lines taken from
*                   Rev. Part. Prop. 1994: Phys. Rev. D50, Pt. I, 1173.
*                 Branching ratios between particles in one line
*                 calculated from CLEBSCH-GORDAN coefficients
*                 N* can be:
*                 name     PDG code      name     PDG code
*                 N(1440)+    12212  N(1440)0    12112
*                 N(1520)+     2124  N(1520)0     1214
*                 N(1535)+    22212  N(1535)0    22112
*                 N(1650)+    32212  N(1650)0    32112
*                 N(1675)+     2216  N(1675)0     2116
*                 N(1680)+    12216  N(1680)0    12116
*                 N(1700)+    22124  N(1700)0    21214
*                 N(1710)+    42212  N(1710)0    42112
*                 N(1720)+    32124  N(1720)0    31214
*
*     Input:      INSTAR:  Number of N* in /PYJETS/
*     ------      
*
*     Output:     
*     -------     
*
*     Called by:  PRODIFFVM
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
*               DECRHO       rho decay
*               DECDEL       Delta decay
*               DECETA       eta decay
*
************************************************************************
#include "rglujets.inc"
#include "rgludat2.inc"
#include "rgpdiss.inc"
#include "rgherup.inc"
****************************************************************
      Double Precision xm1,xw1

      INTEGER INSTAR, IDEC, NCODES,NDECAY,I,J,NNSTAR,ICHECK
      INTEGER NCHARG,NDEC1,NDEC2,NSUM,ITER,IPART,IOK,IHELI,IDUMMY
      REAL R,RANBW
	Double Precision draprn
      EXTERNAL draprn

      DATA ICHECK/0/
      PARAMETER (NCODES = 9)
      PARAMETER (NDECAY = 6)
      INTEGER M
      DOUBLE PRECISION beta(3)

* Legal N*+ codes
      INTEGER NPLCOD (1:NCODES)
     +        /12212, 2124, 22212, 32212, 2216, 12216, 22124,
     +         42212, 32124/
* Legal N*0 codes
      INTEGER N0COD (1:NCODES)
     +        /12112, 1214, 22112, 32112, 2116, 12116, 21214,
     +         42112, 31214/
* Cummulative branching ratios for the channels
* N* -> N pi, Delta pi, N rho, N (pi pi)_I=0, N eta, Lambda K
* for the different mass states
      REAL DECCHN (1:NDECAY-1, 1:NCODES)
     +        /0.67, 0.92, 0.92, 1.00, 1.00,
     +         0.60, 0.80, 1.00, 1.00, 1.00,
     +         0.50, 0.50, 0.50, 0.50, 1.00,
     +         0.75, 0.80, 0.90, 0.90, 0.90,
     +         0.45, 1.00, 1.00, 1.00, 1.00,
     +         0.67, 0.78, 0.89, 1.00, 1.00,
     +         0.10, 0.45, 0.65, 1.00, 1.00,
     +         0.15, 0.40, 0.55, 0.85, 0.85,
     +         0.15, 0.15, 0.90, 0.90, 0.90/
* Cummulative branching ratios for different particle charges
      REAL CGCOEF (1:2, 1:NDECAY)
     +        /0.33, 1.00,
     +         0.50, 0.83,
     +         0.33, 1.00,
     +         0.67, 1.00,
     +         1.00, 1.00,
     +         1.00, 1.00/
* Decay particles for N*+ decays
      INTEGER  NPLDC1 (1:3, 1:NDECAY), NPLDC2 (1:3, 1:NDECAY),
     +         NPLDC3 (1:3, 1:NDECAY)
      DATA ((NPLDC1(I,J), NPLDC2(I,J), NPLDC3(I,J), I=1, 3), J=1,NDECAY)
* N* -> N pi
     +        /2212,  111,    0,
     +         2112,  211,    0,
     +            0,    0,    0,
* N* -> Delta pi
     +         2224, -211,    0,
     +         2214,  111,    0,
     +         2114,  211,    0,
* N* -> N rho
     +         2212,  113,    0,
     +         2112,  213,    0,
     +            0,    0,    0,
* N** -> N (pi pi)_I=0
     +         2212,  211, -211,
     +         2212,  111,  111,
     +            0,    0,    0,
* N* -> N eta
     +         2212,  221,    0,
     +            0,    0,    0,
     +            0,    0,    0,
* N* -> Lambda K
     +         3122,  321,    0,
     +            0,    0,    0,
     +            0,    0,    0/

* Decay particles for N*0 decays
      INTEGER N0DC1 (1:3, 1:NDECAY), N0DC2 (1:3, 1:NDECAY),
     +        N0DC3 (1:3, 1:NDECAY)
      DATA ((N0DC1(I,J), N0DC2(I,J), N0DC3(I,J), I=1, 3), J=1, NDECAY)
* N* -> N pi
     +        /2112,  111,    0,
     +         2212, -211,    0,
     +            0,    0,    0,
* N* -> Delta pi
     +         1114,  211,    0,
     +         2114,  111,    0,
     +         2214, -211,    0,
* N* -> N rho
     +         2112,  113,    0,
     +         2212, -213,    0,
     +            0,    0,    0,
* N** -> N (pi pi)_I=0
     +         2112,  211, -211,
     +         2112,  111,  111,
     +            0,    0,    0,
* N* -> N eta
     +         2112,  221,    0,
     +            0,    0,    0,
     +            0,    0,    0,
* N* -> Lambda K
     +         3122,  311,    0,
     +            0,    0,    0,
     +            0,    0,    0/



      REAL*8 XM (100), PMOM (4, 100), WT, SUM, WIDTH, XMMIN
      REAL*8 PDPCMS (1:4), PDPRS (1:4)

      IDEC  = 0

* if particle already decayed --> error 

      IF (K(INSTAR,1) .NE. 1) THEN
        PRINT *, 'DECNST error: ISTAT (', INSTAR, ') = ', K (INSTAR,1)
        RETURN
      END IF

*
* Look up particle code
*
      NNSTAR = 0
      NCHARG = 0
      DO 1 I = 1, NCODES

*** check if part type is N*+ or N*0

        IF (ABS (K(INSTAR,2)) .EQ. NPLCOD (I)) THEN
          NNSTAR = I
          NCHARG = 1
        ELSE IF (ABS (K(INSTAR,2)) .EQ. N0COD (I)) THEN
          NNSTAR = I
          NCHARG = 0
        END IF
 1    CONTINUE

c      PRINT *, 'NNSTAR = ', NNSTAR, ', NCHARG = ', NCHARG
      IF (NNSTAR .EQ. 0) THEN
        PRINT *, '*** DECNST error: Don''t know particle ',
     +           K(INSTAR,2)
c        CALL ERRLOG (201, 'S: DECNST: Wrong particle!')
        RETURN
      END IF

*
* Choose decay channel
*
 7    CONTINUE

      R = draprn ()
      NDEC1 = 1
 2    CONTINUE
      if(NDEC1.GT.5) goto 3
        IF (NDEC1 .EQ. NDECAY .OR. DECCHN (NDEC1, NNSTAR) .GE. R) GOTO 3
        NDEC1 = NDEC1 + 1
      GO TO 2
 3    CONTINUE

c      IF (NDEC1 .LT. NDECAY) THEN
c        PRINT *, 'R = ', R, ', DECCHN (', NDEC1, ', ', NNSTAR, ') = ',
c     +           DECCHN (NDEC1, NNSTAR)
c      ELSE
c        PRINT *, 'R = ', R, ', DECCHN (', NDEC1-1, ', ', NNSTAR, ') = ',
c     +           DECCHN (NDEC1-1, NNSTAR)
c      END IF

      R = draprn ()
      NDEC2 = 1
 4    CONTINUE
      if(NDEC2.gt.2) goto 5
c      write(6,*) ' ndec ', ndec2,ndec1
        IF (NDEC2 .EQ. 3 .OR. CGCOEF (NDEC2, NDEC1) .GE. R) GO TO 5
        NDEC2 = NDEC2 + 1
      GO TO 4
 5    CONTINUE

*
* Fill the mass array XM
*

      IDEC  = N + 1
c	write(6,*) ' decnst IDEC',IDEC,NDEC2,N,NSUM

      IF (NCHARG .EQ. 1) THEN
        K(IDEC,2)   = NPLDC1 (NDEC2, NDEC1)
        K(IDEC+1,2) = NPLDC2 (NDEC2, NDEC1)
        K(IDEC+2,2) = NPLDC3 (NDEC2, NDEC1)
      ELSE
        K(IDEC,2)   = N0DC1 (NDEC2, NDEC1)
        K(IDEC+1,2) = N0DC2 (NDEC2, NDEC1)
        K(IDEC+2,2) = N0DC3 (NDEC2, NDEC1)
      END IF
      IF (K(INSTAR,2) .LT. 0) THEN
        K(IDEC,2) = -K(IDEC,2)
        IF (K(IDEC+1,2) .NE. 111 .AND. K(IDEC+1,2) .NE. 221) THEN
          K(IDEC+2,2) = -K(IDEC+2,2)
          K(IDEC+3,2) = -K(IDEC+3,2)
        END IF
      END IF


c Decay to K0: 50% K0S, 50% K0L
      IF (ABS (K(IDEC+1,2)) .EQ. 311) THEN
        IF (draprn () .LT. 0.50) THEN
          K(IDEC+1,2) = 130
        ELSE
          K(IDEC+1,2) = 310
        END IF
      END IF

      IF (K(IDEC+2,2) .EQ. 0) THEN
        NSUM = 2
c       PRINT *, ' Decay to ', K(IDEC,2), K(IDEC+1,2)
      ELSE
        NSUM = 3
c       PRINT *, ' Decay to ', K(IDEC,2), K(IDEC+1,2), K(IDEC+2,2)
      END IF

      ITER = 0
c	write(6,*) ' decnst: nsum,instar ',nsum,instar,idec
c	call pylist(1)
 9    CONTINUE
      ITER = ITER + 1
      SUM = 0D0
      DO 10  I = 1, NSUM
cccc        IPART = K(instar+I,2)
        IPART = K(idec-1+I,2)
	  call pdismass(IPART,XM1,XW1)
        WIDTH = XW1
        XM (I) = XM1
c	  write(6,*) ' decnst ', width,xm(i)
c        write(*,*) XM(i)
        XMMIN = XM (I) - 2D0*WIDTH
        IF (ABS (IPART) .EQ. 1114 .OR. ABS (IPART) .EQ. 2114 .OR.
     +      ABS (IPART) .EQ. 2214 .OR. ABS (IPART) .EQ. 2224) THEN
          XMMIN = 1.1D0
        END IF
        IF (WIDTH .GT. 1D-6) THEN
          XM (I) = RANBW (XM (I), WIDTH, XMMIN, XM (I) + 2D0*WIDTH)
        END IF
c	  write(6,*) ' sum,xm ',sum,xm(i),width
        SUM = SUM + XM (I)
 10   CONTINUE


      IF (SUM .GT. P(INSTAR,5) .AND. ITER .LT. 100) GO TO 9
      IF (SUM .GT. P(INSTAR,5)) GO TO 7

*
* Call rgrambo: distribute particles uniformely in phase space
*
         
       do 15 i=1,3
           beta(i)=P(instar,i)/P(instar,4)       
 15     continue

       CALL PYROBO(instar,instar,0.d0,0.d0,-beta(1),-beta(2) 
     &                                ,-beta(3))


       CALL rgrambo (NSUM, P(INSTAR,5), XM, PMOM, WT, 1, IOK)
       IF (IOK .NE. 0) THEN
	    write(6,*) ' decnst: RGRAMBO fails '
c         CALL ERRLOG (202, 'S: DECNST: rgrambo failed!')
       END IF

* Update /PYJETS/
*
 
       do 20 i=1,nsum            
          P(N+i,5)=0
          do 40 j=1,4      
             P(N+i,j)=PMOM(j,i)
c             P(N+i,5)=P(N+i,5)+PMOM(j,i)
             P(N+i,5)=XM(i)
 40       continue
c          P(N+i,4)=PMOM(4,i)
c          P(N+i,5)=P(N+i,4)-P(N+i,5)

 20    continue   
 
c       write(*,*) ' decnst: p(instar) ', P(instar,1),P(instar,2),
c     &              P(instar,3),P(instar,4),P(instar,5)
c       do i=1,nsum
c       write(*,*) ' decnst: p(n+i) ',i,P(N+i,1),P(N+i,2),P(N+i,3),
c     & P(N+i,4),P(N+i,5) 
c       enddo

       CALL PYROBO(instar,instar,0.d0,0.d0,beta(1),beta(2) 
     &                                ,beta(3))
       

       M=N

      K(INSTAR,1) = 11
      K(INSTAR,4) = N + 1
      K(INSTAR,5) = N + NSUM
      M=N
      N=N+NSUM
c	write(6,*) ' decnst ',instar,K(INSTAR,4),K(INSTAR,5),N,NSUM
      DO 50, I = 1, NSUM
        P (N + I,5) = XM (I)
        K (M + I,1) = 1
        K (M + I,3) = INSTAR
        K (M + I,4) = 0
c        K (M + I,5) = 0
 50   CONTINUE

*       N = N + NSUM

       CALL PYROBO(M+1,N,0.d0,0.d0,beta(1),beta(2),beta(3))
      
c       write(*,*) 'rgrambo dcnst'
c
c       call pylist(1)

*
* Perform Delta, pi0, rho, eta decay
*

      IF (ABS (K(IDEC,2)) .EQ. 1114 .OR.
     +    ABS (K(IDEC,2)) .EQ. 2114 .OR.
     +    ABS (K(IDEC,2)) .EQ. 2214 .OR.
     +    ABS (K(IDEC,2)) .EQ. 2224) THEN
 	   K(IDEC,1)=1
      if(ijetset.eq.1) CALL DECDEL (IDEC)
      END IF

      IF (K(IDEC+1,2) .EQ. 111) THEN
	   K(IDEC+1,1)=1	
         if(ijetset.eq.1) CALL DECPI0 (IDEC+1)
      ELSE IF (K(IDEC+1,2) .EQ. 113 .OR.
     +         ABS (K(IDEC+1,2)) .EQ. 213) THEN
        IHELI = INT (3.0*draprn()) - 1
	   K(IDEC+1,1)=1
        if(ijetset.eq.1) CALL DECRHO (IDEC+1, IHELI, IDUMMY)
      ELSE IF (K(IDEC+1,2) .EQ. 221) THEN
	   K(IDEC+1,1)=1
        if(ijetset.eq.1) CALL DECETA (IDEC+1)
      END IF

      IF (K(IDEC+2,2) .EQ. 111) THEN
	   K(IDEC+2,1)=1
        if(ijetset.eq.1) CALL DECPI0 (IDEC+2)
      END IF
	
      if(ijetset.eq.0.and.ILHA.eq.0) call pyexec
      
      K(INSTAR,1) = 11
c      call pylist(2)
      RETURN
      END

