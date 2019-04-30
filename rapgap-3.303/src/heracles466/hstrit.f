C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSTRIT(F,X,NDIM,NCONT,XI,NDO)
C
C  AUTHOR      : J. VERMASEREN
C                MODIFIED 2/88 BY HJM
C                LAST CHANGE 23/07/90  BY HJM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION HSTRIT
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      DIMENSION X(NDIM),XI(50,NDIM)
      DIMENSION Z(10),R(20),NCALL(20)
      DATA NCALL /20*0/
      EXTERNAL F
C
      IF(NCALL(NCONT).EQ.0) THEN
        NCALL(NCONT)=1
        R(NCONT)=NDO**NDIM
      ENDIF
C
      W=R(NCONT)
      DO 4 I=1,NDIM
         XX=X(I)*NDO
         J=XX
         JJ=J+1
         Y=XX-J
         IF(J.LE.0)THEN
            DD=XI(1,I)
         ELSE
            DD=XI(JJ,I)-XI(J,I)
         ENDIF
         Z(I)=XI(JJ,I)-DD*(1.-Y)
         W=W*DD
4     CONTINUE
C
      FZ=F(Z)
      HSTRIT=W*FZ
C
      IF(HSTRIT.LE.0D0.AND.IPRINT.GE.5) THEN
        WRITE(LUNOUT,'(A,5(1PE12.4)/15X,5(1PE12.4))')
     *        ' HSTRIT - Z(I):',Z
        WRITE(LUNOUT,'(A,1PE12.4)')  ' HSTRIT -  W  :',W
        WRITE(LUNOUT,'(A,1PE12.4)')  ' HSTRIT - F(Z):',FZ
        WRITE(LUNOUT,'(A/50(10(1PE12.4)/))')  ' HSTRIT - XI  :',XI
      ENDIF
C
      RETURN
      END
