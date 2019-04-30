C==========================================================
      SUBROUTINE HSRNTE(IO)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
C
C     REAL*4 UU(97)
C     REAL*4 U(6),X(6),D(6)
C
      DIMENSION UU(97)
      DIMENSION U(6),X(6),D(6)
      DATA U / 6533892.0D0 , 14220222.0D0 ,  7275067.0D0 ,
     &         6172232.0D0 ,  8354498.0D0 , 10633180.0D0 /
C
      CALL HSRNOU(UU,CC,CCD,CCM,II,JJ)
      CALL HSRNST(12,34,56,78)
      DO 10 II1 = 1,20000
10      XX      = HSRNDM(-1)
      SD        = 0.0D0
      DO 20 II2 = 1,6
        X(II2)  = 4096.D0*(4096.D0*HSRNDM(-1))
        D(II2)  = X(II2)-U(II2)
20      SD      = SD+D(II2)
      CALL HSRNIN(UU,CC,CCD,CCM,II,JJ)
      IF ( IO.EQ.1 .OR. SD.NE.0D0 )
     &   WRITE(LUNTES,50) (U(I),X(I),D(I),I=1,6)
      RETURN
50    FORMAT('  === TEST OF THE RANDOM-GENERATOR ===',/,
     &       '    EXPECTED VALUE    CALCULATED VALUE     DIFFERENCE',/,
     &       6(F17.1,F20.1,F15.3,/),
     &       '  === END OF TEST ;',
     &       '  GENERATOR HAS THE SAME STATUS AS BEFORE CALLING HSRNTE')
C==END OF RANDOM GENERATOR PACKAGE==========================
      END
