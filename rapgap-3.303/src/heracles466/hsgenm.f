C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSGENM(F,NCONT,NDIM,NEVENT,ICONTI,FFMAX,FMAX,GLAST,
     &                  FFMOLD,FMLOLD,DNCGLO,DNCLOC,LLOCAL,LGLOB,
     &                  NTOT,NM,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,
     &                  XI,NDO,MBIN,NREG)
C
C  ORIGINAL AUTHOR      : J. VERMASEREN
C                MODIFIED 2/88 HJM
C                LAST CHANGE   7/12/89  HJM
C                LAST CHANGE  19/03/90  HJM
C                MAJOR MODIFICATION  12/04/90  HJM
C                LAST CHANGE  23/07/90  HJM
C                MODIFIED 12/02/97 HS
C
C---------------------------------------------------------------------
C
C   VARIABLES:
C   **********
C
C   F :          ACTUAL DISTRIBUTION FUNCTION FOR SAMPLING
C
C   NCONT :      ACTUAL NUMBER OF THE CORRESPONDING CONTRIBUTION
C                TO THE MATRIX ELEMENT
C
C   NDIM  :      DIMENSION OF THE COORDINATE VECTOR
C
C   NEVENT:      NUMBER OF EVENTS TO BE SAMPLED
C
C   ICONTI = 1   NEW SEQUENCE OF EVENTS TO BE SAMPLED,
C                I.E. INFORMATION FROM POTENTIAL PREVIOUS SAMPLING RUNS
C                     IS TO BE DISCARDED
C          = 2   EXTENSION OF AN EXISTING EVENT SAMPLE
C
C   FFMAX :      ESTIMATED GLOBAL MAXIMUM
C
C   FMAX  :      ESTIMATES FOR LOCAL MAXIMA (FIELD)
C
C   GLAST :      LAST CALCULATED FUNCTION VALUE
C                (NEEDED FOR CONTINUED SAMPLING, ICONTI=1)
C
C   NTOT  :      EFFECTIVE TOTAL NUMBER OF TRIALS FOR SAMPLING
C
C   NM(I) :      NUMBER OF TRIALS IN REGION I
C
C   NCALL :      TOTAL NUMBER OF ACTUAL FUNCTION EVALUATIONS
C
C   MCALL1 :     NUMBER OF FUNCTION EVALUATIONS DURING CORRECTION
C                PROCEDURE FOR WRONG LOCAL MAXIMUM
C
C   MCALL2 :     NUMBER OF FUNCTION EVALUATIONS DURING CORRECTION
C                PROCEDURE FOR WRONG GLOBAL MAXIMUM
C
C   IBIMAX :     NUMBER OF THE BIN CONTAINING THE CURRENT GLOBAL MAXIMUM
C
C   FFMOLD :     GLOBAL MAXIMUM FROM LAST SAMPLING RUN
C                (NOT YET CORRECTED EVEN IF NECESSARY!)
C
C   FMLOLD :     LOCAL MAXIMUM FOR THE BIN SELECTED FOR THE LAST EVENT
C                IN THE PREVIOUS SAMPLING RUN
C                (NOT YET CORRECTED EVEN IF NECESSARY!)
C
C   DNCGLO :     NUMBER OF EVENTS TO BE ADDITIONALLY SAMPLED
C                TO CONTINUE THE CORRECTION PROCEDURE
C                FOR WRONG ESTIMATION OF GLOBAL MAXIMUM
C                (INFORMATION FROM POTENTIAL PREVIOUS RUN)
C
C   DNCLOC :     NUMBER OF EVENTS TO BE ADDITIONALLY SAMPLED
C                IN REGION **JCOR** TO CONTINUE THE CORRECTION PROCEDURE
C                FOR WRONG ESTIMATION OF LOCAL MAXIMUM
C                (INFORMATION FROM POTENTIAL PREVIOUS RUN)
C
C   LGLOB :      (LOGICAL) IF TRUE CORRECTION FOR WRONG ESTIMATE OF
C                GLOBAL MAXIMUM NECESSARY
C
C   LLOCAL :     (LOGICAL) IF TRUE CORRECTION FOR WRONG ESTIMATE OF
C                LOCAL MAXIMUM NECESSARY
C
C   XI    :      ACTUAL SUBDIVISION OF NDIM-DIMENSIONAL COORDINATE SPACE
C                AS PROVIDED BY VEGAS
C
C   NDO   :      NUMBER OF INTERVALS PER COORDINATE AXIS FROM VEGAS
C                FROM VEGAS-SUBDIVISION
C
C   MBIN  :      NUMBER OF BINS PER AXIS FOR SUBVOLUMES USED IN EVENT
C                SAMPLING TO INCREASE THE SAMPLING EFFICIENCY
C
C---------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      PARAMETER (NDIMX=10)
      DIMENSION X(NDIMX),N(NDIMX),FMAX(NREG),NM(NREG),XI(NDO,NDIM)
      real*8 ntot,ntold
      LOGICAL LGLOB,LLOCAL
C----------------------------------------------------------------------
C
C  COUNTERS:
C  *********
C
C  NTOT  :   CURRENT TOTAL NUMBER OF TRIALS INCLUDING THE ONES FROM
C            EARLIER RUNS FOR CONTINUED SAMPLING
C  NN    :   CURRENT NUMBER OF TRIALS IN THE ACTUAL RUN
C  NCALL :   TOTAL NUMBER OF ACTUAL FUNCTION CALLS
C  MCALL1:   TOTAL NUMBER OF FUNCTION CALLS DURING CORRECTION PROCEDURE
C            FOR UNDERESTIMATED LOCAL MAXIMA
C  MCALL2:   TOTAL NUMBER OF FUNCTION CALLS DURING CORRECTION PROCEDURE
C            FOR UNDERESTIMATED GLOBAL MAXIMUM
C  NEV   :   CURRENT NUMBER OF ACCEPTED EVENTS IN THE ACTUAL RUN
C  MEV1  :   CURRENT NUMBER OF EVENTS ACCEPTED IN THE ACTUAL RUN
C            DURING CORRECTION PROCEDURE FOR UNDERESTIMATED LOCAL MAXIMA
C  MEV2  :   CURRENT NUMBER OF EVENTS ACCEPTED IN THE ACTUAL RUN
C            DURING CORRECTION PROCEDURE FOR UNDERESTIMATED GLOBAL MAX.
C
C-----------------------------------------------------------------------

      IF(IPRINT.GE.2) THEN
        WRITE(LUNTES,'(//A/A,5X,2I3,I8,I3,2(1PE13.6))')
     +     ' ENTRY HSGENM ',' NCONT,NDIM,NEVENT,ICONTI,FFMAX, GLAST',
     +                        NCONT,NDIM,NEVENT,ICONTI,FFMAX,GLAST
        WRITE(LUNTES,'(5X,A/5X,4I8,4I5)')
     &     ' NTOT,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN,NREG :',
     &       NTOT,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN,NREG
        WRITE(LUNTES,'(5X,A/5X,4(1PD12.5),2L4)')
     &     ' FFMOLD, FMLOLD, DNCGLO, DNCLOC, LLOCAL, LGLOB ',
     &       FFMOLD,FMLOLD,DNCGLO,DNCLOC,LLOCAL,LGLOB
        WRITE(LUNTES,'(A/50(10(1PE12.4)/))')  ' HSGENM - XI  :',XI
      ENDIF

C...Basic initialization
      NN=0
      NEV=0
      MEV1=0
      MEV2=0
      AMI=1D0/FLOAT(MBIN)

      IF(ICONTI.EQ.1) THEN
C...New start of event sampling discarding potentially existing 
C...information on trial numbers and corrections from previous runs
        LGLOB=.FALSE.
        LLOCAL=.FALSE.
        DNCGLO=0D0
        DNCLOC=0D0
        JCOR=0
        NTOT=0
        NCALL=0
        MCALL1=0
        MCALL2=0
        DO 1 I=1,NREG
          NM(I)=0
    1   CONTINUE

      ELSEIF(ICONTI.EQ.2) THEN
C...Restart for continued sampling
        G=GLAST
        GLOMAX=FFMAX
C...Reset parameter for random evaluation of coordinate vector, 
C...needed in case of continued correction
         JJ=JCOR - 1
         DO 3 K=1,NDIM
            JJJ=JJ/MBIN
            N(K)=JJ-JJJ*MBIN
            JJ=JJJ
  3      CONTINUE
C...Test whether correction of local maximum has to be continued
        IF(LLOCAL) GOTO 290
C...Test whether correction of global maximum has to be continued
        IF(LGLOB) GOTO 390
C...Test whether last function value exceeds local/global maxima
        J=JCOR
        GOTO 190

      ELSE
        WRITE(LUNOUT,'(/A,I5/A,I5/A)')
     &       ' *** HSGENM: WRONG INPUT FOR ICONTI, ICONTI=',ICONTI,
     &       ' ***         CALLED WITH NCONT=',NCONT,
     &       ' *** EXECUTION STOPPED'
        stop
      ENDIF

  100 CONTINUE
C...Entry after correctino procedures for local/global maxima
      LLOCAL=.FALSE.
      LGLOB=.FALSE.
      DNCGLO=0D0
      DNCLOC=0D0

  110 CONTINUE
C...Unrestricted event sampling
      NTOT=NTOT + 1
      NN=NN + 1
      J=HSRNDM(-1)*FLOAT(NREG) + 1D0
      Y=HSRNDM(-1)*FFMAX
      NM(J)=NM(J)+1
C...Rejection on the basis of estimated local maxima
      IF(Y.GT.FMAX(J)) GOTO 110
      JJ=J-1
      DO 10 K=1,NDIM
        JJJ=JJ/MBIN
        N(K)=JJ-JJJ*MBIN
        X(K)=(HSRNDM(-1)+N(K))*AMI
        JJ=JJJ
 10   CONTINUE
      G=HSTRIT(F,X,NDIM,NCONT,XI,NDO)
      NCALL=NCALL + 1
C...Rejection on the basis of actual function value
      IF (Y.GE.G) GOTO 110

C...Event accepted
      CALL HSACPT(NCONT)
      NEV=NEV + 1
      IF(NEV.GE.NEVENT) THEN
C...Prepare potential correction for further runs if the last
C...calculated function value is larger than the local/global maximum
        JCOR=J
        GOTO 400
      ENDIF

C*********************************
C...Entry for restart with continued sampling. No further correction
C...from previous runs necessary: test the last function value
C...evaluated in the prevsious runs.
 190  CONTINUE
      IF (G.LE.FMAX(J)) GOTO 110
      JCOR=J

      IF (JCOR.EQ.IBIMAX) THEN
C...Correction of global maximum in bin IBIMAX
        LGLOB=.TRUE.
        GLOMAX=G
        FFMOLD=FFMAX
        GOTO 300
      ENDIF

C...Correct for underestimation of local maximum in bin J.
C...Correction for global maximum later at 300
      LLOCAL=.TRUE.
      FMLOLD=FMAX(J)
      GLOMAX=FFMAX

 200  CONTINUE
C...Entry for correction loop with continued correction of maximum
      IF(G.GT.GLOMAX) THEN
C...New local maximum larger than current global one
        GLOMAX=G
        IF(.NOT.LGLOB) THEN
C...LGLOB=.TRUE. possible at this point for repeated correction in 
C...the same bin
          LGLOB=.TRUE.
          FFMOLD=FFMAX
          DNCLOC=DNCLOC + (NM(JCOR)-1)*(FFMOLD-FMAX(JCOR))/FFMOLD
        ENDIF
      ELSE
C...New value smaller than global maximum
        DNCLOC=DNCLOC + (NM(JCOR)-1)*(G-FMAX(JCOR))/FFMAX
      ENDIF

 201  CONTINUE
      IF(IPRINT.GE.4) THEN
        WRITE(LUNTES,'(/A,I2/A,I5,A,2(1PD14.5))')
     &       ' HSGENM / NCONT=', NCONT,
     &       ' CORRECTION OF LOCAL MAXIMUM IN BIN', JCOR,
     &       ' / OLD-NEW :',FMAX(JCOR),G
      ENDIF
      FMAX(JCOR)=G

 210  CONTINUE
C...Entry for correction loop without continued correction of the 
C...local maximum
      IF(DNCLOC.LT.1D0) THEN
        IF(HSRNDM(-1).GT.DNCLOC) THEN
          DNCLOC=0D0
          LLOCAL=.FALSE.
          GOTO 300
        ENDIF
        DNCLOC=1D0
      ENDIF
      DNCLOC=DNCLOC - 1D0
      DO 211 K=1,NDIM
        X(K)=(HSRNDM(-1)+N(K))*AMI
 211  CONTINUE
      G=HSTRIT(F,X,NDIM,NCONT,XI,NDO)
      MCALL1=MCALL1 + 1
      NCALL=NCALL+1
C...Events accepted only if function value larger than old maximum, 
C...i.e., falling into the bin where events might have been discarded
C...because of wrong maximum
      IF (G.LT.FMLOLD) GOTO 210
      IF (G.LT.(HSRNDM(-1)*FMAX(JCOR))) GOTO 290
      CALL HSACPT(NCONT)
      MEV1=MEV1 + 1
      NEV=NEV + 1
      IF(NEV.GE.NEVENT) THEN
C...Remember potentially outstanding correction for wrong estimation
C...of global maximum
        IF (LGLOB) THEN
          DNCGLO=DNCGLO + NTOT*(GLOMAX-FFMAX)/FFMOLD
          FFMAX=GLOMAX
          IF(DNCGLO.LE.0D0) THEN
            WRITE(LUNTES,'(2A/A,3(1PE13.5))')
     &         ' HSGENM - ERROR AT EXIT: WRONG DEFINITION OF DNCGLO',
     &         ' (BEFORE 290 CONTINUE)', ' G, GLOMAX, FFMAX',
     &          G, GLOMAX, FFMAX
            WRITE(LUNTES,'(A,4I5,2L4/2(A,5(1PE20.12)/))')
     &       ' JCOR, J, IBIMAX, ICONTI, LGLOB, LLOCAL',
     &         JCOR, J, IBIMAX, ICONTI, LGLOB, LLOCAL,
     &       ' FFMAX, GLOMAX, FFMOLD, G, GLAST',
     &         FFMAX, GLOMAX, FFMOLD, G,GLAST,
     &       ' FMAX(JCOR), DNCGLO, DNCLOC, FMLOLD',
     &         FMAX(JCOR), DNCGLO, DNCLOC, FMLOLD
            WRITE(LUNTES,'(5X,A/5X,I20,3I10,4I5)')
     &       ' NTOT,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN,NREG :',
     &         NTOT,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN,NREG
          ENDIF
        ENDIF
        IF (DNCLOC.LE.0D0) LLOCAL=.FALSE.
        GOTO 400
      ENDIF

C****************************
C...Entry for restart with continued sampling / correction for wrong
C...local maximum to be continued: test the function value from the 
C...previous run
 290  CONTINUE
      IF (G.LE.FMAX(JCOR)) THEN
        GOTO 210
      ELSE
C...Further correction of estimated local maximum. New local maximum
C...could become global one
        GOTO 200
      ENDIF

 300  CONTINUE
      IF (.NOT.LGLOB) GOTO 100
C...Correction for underestimation of global maximum
      IBIMAX=JCOR
      NTOLD=NTOT

 320  CONTINUE
C...Entry if repeated correction of global maximum is necessary
      DNCGLO=DNCGLO + NTOLD*(GLOMAX-FFMAX)/FFMOLD
      IF(IPRINT.GE.3) THEN
        WRITE(LUNTES,'(/A,I2/A,I5,A,2(1PD14.5))')
     &       ' HSGENM / NCONT=', NCONT,
     &       ' CORRECTION OF GLOBAL MAXIMUM IN BIN', JCOR,
     &       ' / OLD-NEW :',FFMAX,GLOMAX
      ENDIF
      IF(GLOMAX.LT.FFMAX.OR.DNCGLO.LE.0D0) THEN
        WRITE(LUNOUT,'(/A,I2,A/A,4I5,2L4/2(A,5(1PE20.12)/))')
     &       ' HSGENM / NCONT=', NCONT,
     &       ' : WRONG CORRECTION OF GLOBAL MAXIMUM',
     &       ' JCOR, J, IBIMAX, ICONTI, LGLOB, LLOCAL',
     &         JCOR, J, IBIMAX, ICONTI, LGLOB, LLOCAL,
     &       ' FFMAX, GLOMAX, FFMOLD, G, GLAST',
     &         FFMAX, GLOMAX, FFMOLD, G,GLAST,
     &       ' FMAX(JCOR), DNCGLO, DNCLOC, FMLOLD',
     &         FMAX(JCOR), DNCGLO, DNCLOC, FMLOLD
            WRITE(LUNOUT,'(5X,A/5X,2I22,3I10,4I5)')
     &       ' NTOT,ntold,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN :',
     &         NTOT,ntold,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN
            WRITE(LUNOUT,'(/A)')
     &       ' **** EXECUTION STOPPED '
            STOP
      ENDIF
      FFMAX=GLOMAX
      FMAX(JCOR)=GLOMAX

 310  CONTINUE
      IF(DNCGLO.LT.1D0) THEN
        IF(HSRNDM(-1).GT.DNCGLO) THEN
          DNCGLO=0D0
          LGLOB=.FALSE.
          GOTO 100
        ENDIF
        DNCGLO=1D0
      ENDIF
      DNCGLO=DNCGLO - 1D0
      J=HSRNDM(-1)*FLOAT(NREG) + 1D0
      NM(J)=NM(J)+1
      NTOT=NTOT + 1
      NN=NN + 1
C...Accept additional events only in the bin containing the 
C...underestimated global maximum
      IF (J.NE.JCOR) GOTO 310
      DO 311 K=1,NDIM
        X(K)=(HSRNDM(-1)+N(K))*AMI
 311  CONTINUE
      G=HSTRIT(F,X,NDIM,NCONT,XI,NDO)
      MCALL2=MCALL2 + 1
      NCALL=NCALL + 1

C...Event accepted only if function value larger than old maximum to 
C...correct for too high relative normalization of all other bins.
      IF(G.LE.FFMOLD) GOTO 310
C...Standard rejection procedure (step 2) for events in the required
C...region
      IF (G.GE.(HSRNDM(-1)*FFMAX)) THEN
        CALL HSACPT(NCONT)
        MEV2=MEV2 + 1
        NEV=NEV + 1
        IF(NEV.GE.NEVENT) THEN
          IF (DNCGLO.LE.0D0) LGLOB=.FALSE.
C...Remember potentially outstanding correction for wrong estimation 
C...of the global maximum
          GOTO 400
        ENDIF
      ENDIF

C****************************
C...Entry for restart with continued sampling / correction of wrong
C...global maximum to be continued: Test the function value of the 
C...previous run
 390  CONTINUE
      IF (G.GT.FFMAX) THEN
        GLOMAX=G
        GOTO 320
      ELSE
        GOTO 310
      ENDIF

 400  CONTINUE
      GLAST=G

      IF(IPRINT.GE.3) THEN
        WRITE(LUNTES,
     &         '(/A/5X,A/5X,I3,4(1PD12.5)/5X,A/5X,2(1PD12.5),2L4)')
     &  ' *** TEST PRINT HSGENM / EXIT',
     &     ' NCONT, FFMAX, FFMOLD, FMLOLD, GLAST ',
     &       NCONT,FFMAX, FFMOLD,FMLOLD, GLAST,
     &     ' DNCGLO, DNCLOC, LLOCAL, LGLOB ',
     &       DNCGLO, DNCLOC, LLOCAL, LGLOB
        WRITE(LUNTES,'(5X,A/5X,6I8,2I5)')
     &     ' NTOT, NCALL, MCALL1, MCALL2, MEV1, MEV2, IBIMAX, JCOR',
     &       NTOT,NCALL,MCALL1,MCALL2,MEV1,MEV2,IBIMAX,JCOR
        WRITE(LUNTES,'(/A)') ' *****  ARRAY NM(IBIN)'
        WRITE(LUNTES,'(10I8)')  (NM(I),I=1,NREG)
        WRITE(LUNTES,'(/A)') ' *****  LOCAL MAXIMA FMAX(I) '
        WRITE(LUNTES,'(3(I5,1PE15.6))')  (I,FMAX(I),I=1,NREG)
      ENDIF
      RETURN
      END
