      FUNCTION FXN1(X,WEIGHT)
      IMPLICIT NONE
      Double Precision X,WEIGHT,WEIG,FXN1
      DIMENSION X(20)
      Integer Ncall,I,J
#include "rglujets.inc"
#include "rgparas.inc"
#include "rgpart.inc"
#include "rgdisdif.inc"
      DATA NCALL/0/
      DO 10  I=1,20
         DO 10 J=1,5
            K(I,J) = 0
   10 P(I,J) = 0.0D0
      N=0
      WEIG=0.D0
      FXN1=0.D0
      IF(IDIR.EQ.0) THEN
         IF(IPRO.EQ.10.OR.IPRO.EQ.13.OR.IPRO.EQ.15.OR.IPRO.EQ.16) THEN
            CALL DIFFR1(X,WEIG)
         ELSEIF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
            CALL DIFFR2(X,WEIG)
         ELSEIF(IPRO.EQ.12) THEN
            CALL DIFFR3(X,WEIG)
         ELSEIF(IPRO.EQ.18) THEN
            CALL DIFFR8(X,WEIG)
         ELSEIF(IPRO.EQ.20) THEN
            CALL DIFFR5(X,WEIG)
         ELSEIF(IPRO.EQ.21) THEN
            CALL DIFFR6(X,WEIG)
         ELSEIF(IPRO.EQ.22) THEN
            CALL DIFFR7(X,WEIG)
         ELSEIF(IPRO.EQ.30) THEN
            CALL DIFFR30(X,WEIG)
         ELSEIF(IPRO.EQ.100) THEN
            CALL DIFFR4(X,WEIG)
         ELSEIF(IPRO.EQ.500) THEN
            CALL DIFXV1(X,WEIG)
         ELSEIF(IPRO.EQ.501) THEN
            CALL DIFXVq(X,WEIG)
         ELSE
            WRITE(6,*) ' FXN1: wrong subprocess selected; IPRO = ',
     +      IPRO
            WRITE(6,*) ' FXN1: PROGRM STOPPED '
            STOP
         ENDIF
      ELSEIF(IDIR.EQ.1) THEN
         IF(IPRO.EQ.10.OR.IPRO.EQ.13.OR.IPRO.EQ.15.OR.IPRO.EQ.16) THEN
            CALL DIR1(X,WEIG)
         ELSEIF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
            CALL DIR2(X,WEIG)
         ELSEIF(IPRO.EQ.12) THEN
            CALL DIR3(X,WEIG)
         ELSEIF(IPRO.EQ.18) THEN
            CALL DIR4(X,WEIG)
         ELSEIF(IPRO.EQ.99) THEN
            CALL DIR99(X,WEIG)
         ELSEIF(IPRO.EQ.102) THEN
            CALL DIRh0(X,WEIG)
         ELSEIF(IPRO.EQ.500) THEN
            CALL DIRXV1(X,WEIG)
         ELSEIF(IPRO.EQ.501) THEN
            CALL DIRXVq(X,WEIG)
c         elseif(ipro.ge.503.and.ipro.le.509) then
c           Call DIRVQQ(X,Weig)
           ELSE
           WRITE(6,*) ' FXN1: wrong subprocess selected; IPRO = ',
     +      IPRO
           WRITE(6,*) ' FXN1: PROGRM STOPPED '
           STOP
         ENDIF
      ENDIF
C BOOST BACK TO OVERALL CMS
c      write(6,*) '  FXN1 ',weig,X
      IF(WEIG.EQ.0.0) RETURN
      IF(IPRO.EQ.20.OR.IPRO.EQ.21.OR.IPRO.EQ.22.OR.IPRO.EQ.30) THEN
      ELSE
         IF(IPRO.NE.12.OR.IPRO.EQ.100) CALL PYROBO(0,N,0.D0,0.D0,
     +   DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4), DBCMS(3)/DBCMS(4))
      ENDIF
      ncall = ncall + 1
      IF(WEIG.LE.0.0) then
        write(6,*) ' PROBLEM: FXN1 < 0 ; set to 0 '
        WEIG=0.0
      ENDIF
      FXN1=WEIG
      WEIGHT = WEIG
      RETURN
      END
