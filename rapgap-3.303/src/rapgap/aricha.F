
      SUBROUTINE ARICHA
*#**********************************************************************
*#
*#    SUBROUTINE ARICHA
*#
*# PURPOSE: interpret the BOS text bank GARI to change the default
*#          parameters of ARIADNE.
*#
*#
*#**********************************************************************
      character *4   param 
*
	Integer nloop,nmax
	Parameter (nmax=1000)
      character *4   para1 
	Real Rval 
	Integer Ival,Ld,Le
      character*50 Cval 
	Common/steering/Nloop,Ld(nmax),Le(nmax),Ival(nmax),
     &	para1(nmax),Rval(nmax),Cval(nmax)

* ARIADNE parameter common
      COMMON /ARDAT1/ PARA(40), MSTA(40)
      SAVE /ARDAT1/

C---Determine the machine accuacy EPSMAC (D= 2.E-06)
      EPSMAC = 0.5
      DO 33 I= 1, 100
      EPSMAC = EPSMAC * 0.5
      IF ((1.0+EPSMAC) .EQ. 1.0)  GO TO 35
   33 CONTINUE
      EPSMAC = 1.0E-6
   35 EPSMAC = 2.0 * EPSMAC
      WRITE(6,330) EPSMAC
 330  FORMAT(1X,/,' COMPUTED MACHINE ACCURACY = ',G10.2,/)

C---Interpret the text containing the parameters to be changed
      If(Nloop.ge.1) Then
	      Do 30 I=1,Nloop

            LND = Ld(I)
            LNE = Le(I)

            PARAM = Para1(I)
            INT = Ival(I)
            REA = Rval(I)

            IF (PARAM.EQ.'MSTA') THEN
                IF( LND .LE. 0 .OR. LND .GT. 40 )     GO TO 999
                IDEF = MSTA(LND)
                IF(ABS(MSTA(LND)-INT).GT.EPSMAC) THEN
                    MSTA(LND) = INT
                    PRINT 1001, LND,IDEF,INT
                ENDIF
            ENDIF
            IF (PARAM.EQ.'PARA') THEN
                IF( LND .LE. 0 .OR. LND .GT. 40 )     GO TO 999
                RDEF = PARA(LND)
                IF(ABS(PARA(LND)-REA).GT.EPSMAC) THEN
                    PARA(LND) = REA
                    PRINT 1002, LND,RDEF,REA
                ENDIF
            ENDIF

30      CONTINUE

      ENDIF
      RETURN
999   CONTINUE
      PRINT 9999, PARAM, LND, LNE
	   
 1001 FORMAT(' aricha:',
     &'Parameter MSTA(',I4,') changed from default ',I12,' TO ',I12)
 1002 FORMAT(' aricha:',
     &'Parameter PARA(',I4,') changed from default ',E13.6,
     &' TO ',E13.6 )
 9999 FORMAT(' aricha:',/,
     &'**************** ERROR: INVALID ARRAY INDEX ************',/,
     &'**************** PARAM =',A4,'  1. INDEX =',I4,
     &'  2. INDEX =',I4 )
      RETURN

      END
