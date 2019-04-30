      SUBROUTINE PTIME(NAME,J,K)
*#*******************************************************************
*#                                                                  *
*#     SUBROUTINE PTIME(NAME,J,K)                                   *
*#                                                                  *
*#     ADDOPTED FROM ROUTINE UFTIME                                 *
*# PURPOSE: HANDLING OF DIFFERENT TIMING OPTIONS                    *
*#                                                                  *
*#                                                                  *
*# INPUT :  NAME : NAME OF TIME INTERVALL                           *
*#          J = 1 BEGIN, J = 2 END                                  *
*#          K > 0 FLAG FOR INCREMENTATION OF DIFF COUNTERS          *
*#          K = 99 TERMINATION, PRINTOUT                            *
*#                                                                  *
*# OUTPUT : VIA COMMON SFDEBU                                       *
*#                                                                  *
*
      PARAMETER (MAXNAM=50)
      CHARACTER*8 NAME
      CHARACTER*8 LNAME(MAXNAM)
      REAL UFTIMX(MAXNAM,4)

      DATA ISTART/1/
      DATA ICALL/0/
      DATA NLOOP/0/
      DATA IERR/0/

C-----
      ICALL=ICALL+1
      IF(ISTART.EQ.1)THEN
         ISTART=0
         N=MAXNAM*4
         CALL VZERO(UFTIMX,N)
         DO 10 I=1,MAXNAM
   10    LNAME(I) = ' '
      ENDIF

      IF(K.EQ.99) GOTO 50

C--------DETERMINE INDEX OF NAME

      DO 20  I=1,NLOOP
         IF(NAME.EQ.LNAME(I)) GOTO 30
   20 CONTINUE
      I=0
   30 CONTINUE

      IF(I.EQ.0) THEN
         IF(NLOOP.GE.MAXNAM) THEN
            WRITE(*,*)' PTIME ERROR: TOO MANY TIMING OPTIONS'
            IERR=1
            RETURN
         END IF
         NLOOP=NLOOP+1
         I=NLOOP
         LNAME(I)=NAME
      END IF

C        WRITE(*,3000)NAME,I,J,K,NLOOP
C3000    FORMAT(' PTIME: NAME,I,J,K,NLOOP= ',A8,2X,I5,2X,I5,2X,I5,2X,I5)

C--------MODULE TIMING

      CALL TIMED(TIME)
C        WRITE(*,*)' TIME=',TIME
      DO 40  ILOOP=1,NLOOP
         UFTIMX(ILOOP,2) = UFTIMX(ILOOP,2) + TIME * UFTIMX(ILOOP,1)
   40 CONTINUE

C--------START OF MODULE TIMING

      IF(J.EQ.1) THEN
         UFTIMX(I,1) = 1.
         UFTIMX(I,3) = UFTIMX(I,3) + 1.

C--------END OF MODULE TIMING

      ELSE
         UFTIMX(I,1) = 0.
         IF(UFTIMX(I,3).eq.1) UFTIMX(I,4) = UFTIMX(I,2)
      END IF

      RETURN

C--------PRINT OUT OF RESULTS
   50 CONTINUE
      IF(NLOOP.LE.0) GOTO 70
* escape if logical unit set to zero
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
     +' *************************************************************'
      WRITE(*,*)
     +' ********************  GENERAL TIMING   **********************'
      WRITE(*,*)
     +' *************************************************************'
      WRITE(*,*)
     +' *********************** IN SECONDS **************************'

C-----TOO MANY LOOPS ARE REQUIRED
      WRITE(*,*)
      IF(IERR.EQ.1) WRITE(*,*)' ERROR ERROR ERROR  TOO MANY LOOPS'

      WRITE(*,*)
     +' ROUTINE        TIME      #CALLS   TIME/1._CALL    TIME/CALL'
      DO 60 I=1,NLOOP
         IF(UFTIMX(I,3).GT.0.5) WRITE(*,10000) LNAME(I), UFTIMX(I,2),
     +   IFIX(UFTIMX(I,3)), UFTIMX(I,4), (UFTIMX(I,2)-UFTIMX(I,4))/
     +   MAX(1.,UFTIMX(I,3)-1.)
   60 CONTINUE

C     WRITE(*,*)
C     WRITE(*,*)' NUMBER OF RECORDS =',IEVENT
C     WRITE(*,*)' NUMBER OF EVENTS  =',KEVENT

C1000 FORMAT(' ',A8,2X,F12.6,2X,I7,2X,F11.6)
10000 FORMAT('  ',A8,2X,F12.6,2X,I7,2X,F12.6,2X,F12.6)
      WRITE(*,*)' *****************************************************'
   70 CONTINUE
      RETURN
      END

