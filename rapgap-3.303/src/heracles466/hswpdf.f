C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   H. SPIESBERGER 22.03.91
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   SUBROUTINE FOR STRUCTURE FUNCTIONS
C   CALCULATED FROM PARTON DISTRIBUTION FUNCTIONS
C   UPDATED: 28.08.92 FOR LOW Q2
C   UPDATED: 13.08.93 FOR LOW FL
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   FOR CONVENTIONS OF THE USE OF PARTON DISTRIBUTION FUNCTIONS
C   AND STRUCTURE FUNCTIONS AT LOW Q2 SEE THE MANUAL
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSWPDF
C---CHECK CONSISTENCY OF PARTON DISTRIBUTION OR STRUCTURE FUNCTION
C   INPUT AND WRITE CHOSEN OPTIONS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSNUCL/ HNA,HNZ
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      REAL*4          PYSTOP,PYSLAM
      COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
      REAL*4          PAR111,PAR112,PARL11,PARL19
      INTEGER         MST111,MST115
C
C---OPTION FOR STRUCTURE FUNCTION INPUT
C---CHECK CONSISTENCY WITH OTHER INPUT DEFINITIONS
      IOPCHQ=0
      IF (IPDFOP.EQ.0.OR.IFLOPT.GE.1) THEN
        DO 3 IL=13,16
        IF (LPAR(IL).NE.0) THEN
          LPAR(IL)=0
          IOPCHQ=1
        ENDIF
    3   CONTINUE
        IF (IPDFOP.EQ.0)
     *  WRITE(LUNOUT,'(10X,A,I4)')
     *        ' STRUCTURE FUNCTION INPUT: IPDFOP = ',IPDFOP
        IF (IOPCHQ.NE.0) THEN
          WRITE(LUNOUT,'(4(10X,A,/))')
     *    ' WARNING: ONLY LEPTONIC CORRECTIONS AND SELF ENRGIES',
     *    ' CAN BE APPLIED FOR STRUCTURE FUNCTION INPUT',
     *    ' OR IF F_L IS INCLUDED ',
     *    ' LPAR(13...16) ARE SET TO 0'
        ENDIF
        IOPCHQ=0
        IF(INT3(4).GT.0 .OR. ISAM3(4).GT.0) THEN
          INT3(4)=0
          ISAM3(4)=0
          IOPCHQ=1
        ENDIF
        IF (IOPCHQ.NE.0) THEN
          WRITE(LUNOUT,'(3(10X,A/))')
     *    ' WARNING: QUARKONIC BREMSSTRAHLUNG CANNOT BE SEPARATED',
     *    ' FOR STRUCTURE FUNCTION INPUT OR IF F_L IS INCLUDED',
     *    ' INT3(4) AND ISAM3(4) ARE SET TO 0'
        ENDIF
      ENDIF

C---WRITE OPTIONS
      IF (ILIB.EQ.1) THEN
        WRITE(LUNOUT,'(/A/A/)')
     *   '           PARTON DISTRIBUTIONS TAKEN FROM PYSTFU   *****'
     *  ,'           VIA LYSTFU FROM LEPTO 6.5                *****'
        IF (ICODE.NE.0) THEN
        WRITE(LUNOUT,'(A,I5)')
     *   '           WITH IDENTIFICATION CODE = ',ICODE
        ELSE
        WRITE(LUNOUT,'(A)')
     *   '           SEE VALUE OF ILQMOD AND THE MANUAL'
        ENDIF
      ELSEIF (ILIB.EQ.2) THEN
        WRITE(LUNOUT,'(/A)')
     *   '           PARTON DISTRIBUTIONS TAKEN FROM PDFLIB   *****'
        WRITE(LUNOUT,'(A,I5)')
     *   '           WITH IDENTIFICATION CODE IVAL = ',ICODE
      ELSEIF (ILIB.EQ.3) THEN
        WRITE(LUNOUT,'(/A/A/)')
     *   ' *****  WARNING: WRONG CODE FOR ILIB: ILIB=3 NOT    *****'
     *  ,'                 ALLOWED (IS OBSOLETE)              *****'
      ENDIF
      WRITE(LUNOUT,'(/A/)')
     *      ' *****  LOW Q2 MODEL FOR STRUCTURE FUNCTIONS:    *****'
      WRITE(LUNOUT,'(A,I2)')
     *   '           ILQMOD = ',ILQMOD
      WRITE(LUNOUT,'(/8(A/))')
     * '           ILQMOD =  0: UNMODIFIED PARTON DISTRIBUTIONS'
     *,'                  =  1: LOW Q2 SUPPRESSED PDF    '
     *,'                  =  2: BRASSE AND STEIN WITH PDF           '
     *,'                  =  3: ALLM(1997) WITH PDF  '
     *,'                  =  4: BADELEK AND KWIECINSKI WITH PDF '
     *,'                  =  5: DONNACHIE AND LANDSHOFF WITH PDF '
     *,'                       (FOR DETAILS, SEE THE MANUAL )'
     *,'                  = 10: STRUCTURE FUNCTIONS FROM USER ROUTINE' 
      IF (ILQMOD.GE.4.AND.ILQMOD.NE.10) THEN
        WRITE(LUNOUT,'(//3(A/))')
     *' ***** WARNING: MAKE SURE THAT PARTON DISTRIBUTION FUNCTIONS'
     *,'                ARE CHOSEN CONSISTENTLY WITH THE LOW Q2 '
     *,'                BEHAVIOUR OF THE F_2 PARAMETRIZATION '
      ENDIF
      IF (ILQMOD.EQ.4.AND.XMIN.LT.1D-5) THEN
        WRITE(LUNOUT,'(//3(A/))')
     *' ***** WRONG COMBINATION OF STRUCTURE FUNCTION CODE WITH '
     *,'       KINEMATIC LIMITS: BADELEK-KWIECINSKI NOT VALID FOR '
     *,'       X-VALUES BELOW 1E-5. EXECUTION STOPPED '
        STOP
      ENDIF

C---FOR CHARGED CURRENT: NO FL, REMOVE B QUARKS
      IF (INT2(2).GT.0
     *   .OR.INT3(7).GT.0.OR.INT3(8).GT.0.OR.INT3(9).GT.0
     *   .OR.ISAM2(2).GT.0
     *   .OR.ISAM3(7).GT.0.OR.ISAM3(8).GT.0.OR.ISAM3(9).GT.0
     *   ) THEN
        IF (IFLOPT.NE.0) THEN
          IFLOPT=0
          WRITE(LUNOUT,'(//A/)')
     *' ***** LONGITUDINAL STRUCTURE FUNCTION NOT INCLUDED FOR CC *****'
        ENDIF
        IF (NPYMAX.GE.5) THEN
          NPYMAX=4
          WRITE(LUNOUT,'(//A/)')
     *' ***** NUMBER OF FLAVORS: NPYMAX RESET TO 4 (NO TOP IN CC) *****'
        ENDIF
      ENDIF

C---OPTIONS FOR LONGITUDINAL STRUCTURE FUNCTION
C---TO BE USED TOGETHER WITH PARTON DISTRIBUTIONS
      IF (IFLOPT.EQ.0) THEN
        WRITE(LUNOUT,'(//A/)')
     *    ' *****  LONGITUDINAL STRUCTURE FUNCTION NOT INCLUDED *****'
      ELSEIF (IFLOPT.GE.1) THEN
        IF (ILQMOD.GE.2) THEN
          IFLOPT=0
          LQCD=0
          LTM=0
          LHT=0
          WRITE(LUNOUT,'(//4(A/))')
     *      ' *****  WARNING: LONGITUDINAL STRUCTURE FUNCTION *****'
     *     ,'           NOT INCLUDED. '
     *     ,'           INCONSISTENT INPUT: IFLOPT > 0 AND '
     *     ,'           ILQMOD > 1, IFLOPT SET TO 0.'
        ELSE
          LQCD=MOD(IFLOPT,10)
          LTM=MOD(IFLOPT/10,10)
          LHT=IFLOPT/100
          IF (ILQMOD.LE.1) IPDFOP=2
          WRITE(LUNOUT,'(//A/,3(A,I3,/))')
     *      ' *****  LONGITUDINAL STRUCTURE FUNCTION INCLUDED *****'
     *     ,'           QCD CONTRIBUTION TO F_L: LQCD = ',LQCD
     *     ,'           TARGET MASS EFFECTS:      LTM = ', LTM
     *     ,'           HIGHER TWIST:             LHT = ', LHT
          WRITE(LUNOUT,'(//A/,2(A,I5,/),2(A,F10.4,/))')
     *      ' *****  DETERMINATION OF ALPHA_S NEEDED FOR F_L: *****'
     *     ,'           ORDER OF ALPHA_S IN ULALPS: MST111 = ',MST111
     *     ,'           TREATMENT OF SINGULARITY:   MST115 = ',MST115
     *     ,'           FIX ALPHA_S VALUE:          PAR111 = ',PAR111
     *     ,'           LAMBDA IN RUNNING ALPHA_S:  PAR112 = ',PAR112
     *     ,'           ACCUARCY IN FL-INTEGRATION: PARL11 = ',PARL11
     *     ,'           PARAMETER FOR HIGHER TWIST: PARL19 = ',PARL19
          WRITE(LUNOUT,'(//A/A/A)')
     *     ' *****  NOTE: LOW Q2 BEHAVIOUR OF F_L IS DETERMINED ****'
     *    ,'              BY THE VALUE OF ILQMOD, '
     *    ,'              SEE MANUAL FOR DETAILS '
          IF (IPDFOP.LT.2) WRITE(LUNOUT,'(//A/A/A)')
     *      ' *****  WARNING: WITH THIS OPTION NO SEPARATION ****'
     *     ,'        OF FLAVORS FOR THE TOTAL CROSS SECTION '
     *     ,'        ---> DJANGO CAN NOT RUN '
        ENDIF
      ENDIF

C---NUCLEAR TARGET
      IF (HNA.NE.1.OR.HNZ.NE.1) THEN
        IF (ILQMOD.LE.1) IPDFOP=2
        WRITE(LUNOUT,'(//A/)')
     *    ' *****  NUCLEAR TARGET  *****'
        WRITE(LUNOUT,'(/2(A,F5.0,/))')
     *    '           A-NUCLEUS = ',HNA
     *   ,'           Z-NUCLEUS = ',HNZ
      ENDIF
      RETURN
      END
