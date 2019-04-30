C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFLAV(W2,IFL,IFLR)

C...Choose flavour of struck quark and the
C...corresponding flavour of the target remnant jet.

      COMMON /HSCUMS/  CQP(12)
      DOUBLE PRECISION CQP,HSRNDM,R

      NFL=6
   20 R=HSRNDM(-1)*CQP(12)
      DO 30 I=1,2*NFL
      IFL=I
      IF(R.LE.CQP(I)) GOTO 40
   30 CONTINUE
   40 CONTINUE
      IF(MOD(IFL,2).EQ.1) THEN
        IFL=(IFL+1)/2
      ELSE
        IFL=-IFL/2
      ENDIF

      IFLR=-IFL

      RETURN
      END
