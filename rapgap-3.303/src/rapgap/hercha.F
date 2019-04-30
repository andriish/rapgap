
      SUBROUTINE HERCHA

*#**********************************************************************
*#
*#    SUBROUTINE HERCHA
*#
*# PURPOSE: Interpret the seering to change the default
*#          Parameters of HERACLES	
*#          note: only those Parameters can be changed which the user
*#
*#**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      character *4   param 
*
	Integer nloop,nmax
	Parameter (nmax=1000)
      character *4   para 
	Real Rval 
	Integer Ival,Ld,Le
      character*50 Cval 
	Common/steering/Nloop,Ld(nmax),Le(nmax),Ival(nmax),
     &	para(nmax),Rval(nmax),Cval(nmax)


C Commons FOR HERACLES.

      PARAMETER ( NCHN2=1, NCHC2=2)
      PARAMETER ( NCHN31=6, NCHN32=7, NCHN33=8, NCHN34=9)
      PARAMETER ( NCHC31=12, NCHC32=13, NCHC33=14)
      PARAMETER(NDIM2=2,NBIN2=50)
      COMMON /HSTCUT/ THEMIN,CTHMIN,CTHCON
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     +                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW/  SW,CW,SW2,CW2
     +              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     +              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMMON /HSISGM/ TCUTQ,TCUTQS
      PARAMETER(NREG2N=2500)
      LOGICAL LGLO2,LLOC2
      COMMON /HSSNC2/ SIG2,SIG2E,T2GGMA,T2GMAX(NREG2N),
     +                XX2(50,2),
     +                FFGO2,DNCG2,FFLO2,DNCL2,GOLD2,
     +                NM2(NREG2N),NDO2,
     +                NTOT2,NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     +                LGLO2,LLOC2

      PARAMETER(NREG2C=2500)
      LOGICAL LGLO2C,LLOC2C
      COMMON /HSSCC2/ SIG2C,SIG2EC,T2GMAC,T2MAXC(NREG2C),
     +                XX2C(50,2),
     +                FFGO2C,DNCG2C,FFLO2C,DNCL2C,GOLD2C,
     +                NM2C(NREG2C),NDO2C,
     +                NTOT2C,NCAL2C,NCA12C,NCA22C,IBIM2C,JCOR2C,
     +                LGLO2C,LLOC2C
      PARAMETER(NDIM31=5,NBIN31=6,NREG31=7776)
      LOGICAL LGLO31,LLOC31
      COMMON /HSSN31/ SIG31,SIG31E,T31GMA,T31MAX(NREG31),
     +                XX31(50,NDIM31),
     +                FFGO31,DNCG31,FFLO31,DNCL31,GOLD31,
     +                SI31,SI2N31,SWGT31,SCHI31,IT31,
     +                NM31(NREG31),NDO31,
     +                NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31,
     +                LGLO31,LLOC31
      PARAMETER(NDIM32=5,NBIN32=6,NREG32=7776)
      LOGICAL LGLO32,LLOC32
      COMMON /HSSN32/ SIG32,SIG32E,T32GMA,T32MAX(NREG32),
     +                XX32(50,NDIM32),
     +                FFGO32,DNCG32,FFLO32,DNCL32,GOLD32,
     +                SI32,SI2N32,SWGT32,SCHI32,IT32,
     +                NM32(NREG32),NDO32,
     +                NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32,
     +                LGLO32,LLOC32
      PARAMETER(NDIM33=5,NBIN33=6,NREG33=7776)
      LOGICAL LGLO33,LLOC33
      COMMON /HSSN33/ SIG33,SIG33E,T33GMA,T33MAX(NREG33),
     +                XX33(50,NDIM33),
     +                FFGO33,DNCG33,FFLO33,DNCL33,GOLD33,
     +                SI33,SI2N33,SWGT33,SCHI33,IT33,
     +                NM33(NREG33),NDO33,
     +                NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33,
     +                LGLO33,LLOC33
      PARAMETER(NDIM34=5,NBIN34=6,NREG34=7776)
      LOGICAL LGLO34,LLOC34
      COMMON /HSSN34/ SIG34,SIG34E,T34GMA,T34MAX(NREG34),
     +                XX34(50,NDIM34),
     +                FFGO34,DNCG34,FFLO34,DNCL34,GOLD34,
     +                SI34,SI2N34,SWGT34,SCHI34,IT34,
     +                NM34(NREG34),NDO34,
     +                NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34,
     +                LGLO34,LLOC34
      PARAMETER(NDM3CC=5)
      PARAMETER(NBN31C=6,NRG31C=7776)
      LOGICAL LGL31C,LLC31C
      COMMON /HSSC31/ SIG31C,SG31EC,T31GMC,T31MXC(NRG31C),
     +                XX31C(50,5),
     +                FFG31C,DNG31C,FFL31C,DNL31C,GLD31C,
     +                SI31C,S2N31C,SWT31C,SCH31C,IT31C,
     +                NM31C(NRG31C),NDO31C,
     +                NTT31C,NCL31C,NC131C,NC231C,IBM31C,JCR31C,
     +                LGL31C,LLC31C
      PARAMETER(NBN32C=6,NRG32C=7776)
      LOGICAL LGL32C,LLC32C
      COMMON /HSSC32/ SIG32C,SG32EC,T32GMC,T32MXC(NRG32C),
     +                XX32C(50,5),
     +                FFG32C,DNG32C,FFL32C,DNL32C,GLD32C,
     +                SI32C,S2N32C,SWT32C,SCH32C,IT32C,
     +                NM32C(NRG32C),NDO32C,
     +                NTT32C,NCL32C,NC132C,NC232C,IBM32C,JCR32C,
     +                LGL32C,LLC32C
      PARAMETER(NBN33C=6,NRG33C=7776)
      LOGICAL LGL33C,LLC33C
      COMMON /HSSC33/ SIG33C,SG33EC,T33GMC,T33MXC(NRG33C),
     +                XX33C(50,5),
     +                FFG33C,DNG33C,FFL33C,DNL33C,GLD33C,
     +                SI33C,S2N33C,SWT33C,SCH33C,IT33C,
     +                NM33C(NRG33C),NDO33C,
     +                NTT33C,NCL33C,NC133C,NC233C,IBM33C,JCR33C,
     +                LGL33C,LLC33C
      CHARACTER*45 CHNAME
      COMMON /HSNAMC/ CHNAME(20)
      COMMON /HSNUME/ SIGTOT,SIGTRR,SIGG(20),SIGGRR(20),NNEVNT,NEVE(20)
      COMMON /HSVGLP/ NPOVEG,NUMINT,NPHYP
      COMMON /HSRDIO/ ISDINP,ISDOUT
      COMMON /VGRES/  S1,S2,S3,S4
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      REAL*4           PYSTOP,PYSLAM
      COMMON /PYSTFUC/ PYSTOP,PYSLAM,NPYMOD,NPYMAX,NPYMIN
      REAL*4 PAR111,PAR112,PARL11,PARL19
      INTEGER MST111,MST115
	COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
      COMMON /HSWGTC/ IWEIGS



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

            PARAM = Para(I)
            INT = Ival(I)
            REA = Rval(I)
		
		
C---------------this is the Heracles part-------------------------------
            IF( PARAM .EQ. 'POLA' ) THEN
               RDEF = REAL(POLARI)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  POLARI = DBLE(REA)
                  PRINT 10400, RDEF, REA
               END IF
            END IF
            IF( PARAM .EQ. 'ICUT' ) THEN
               IDEF = ICUT
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  ICUT = INT
                  PRINT 10700, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'WMIN' ) THEN
               RDEF = REAL(WMIN)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  WMIN = DBLE(REA)
                  PRINT 11200, RDEF, REA
               END IF
            END IF
            IF( PARAM .EQ. 'XMAX' ) THEN
               RDEF = REAL(XMAX)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  XMAX = DBLE(REA)
                  PRINT 11300, RDEF, REA
               END IF
            END IF
            IF( PARAM .EQ. 'XMIN' ) THEN
               RDEF = REAL(XMIN)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  XMIN = DBLE(REA)
                  PRINT 10800, RDEF, REA
               END IF
            END IF
            IF( PARAM .EQ. 'EGMN' ) THEN
               RDEF = REAL(EGMIN)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  EGMIN = DBLE(REA)
                  PRINT 11400, RDEF, REA
                  IF( EGMIN.GT.0D0 ) IOPEGM = 1
               END IF
            END IF
            IF( PARAM .EQ. 'INT2' ) THEN
               IF( LND .LE. 0 .OR. LND .GT. 5 ) GO TO 30
               IDEF = INT2(LND)
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  INT2(LND) = INT
                  PRINT 11500, LND, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'INT3' ) THEN
               IF( LND .LE. 0 .OR. LND .GT. 15 ) GO TO 30
               IDEF = INT3(LND)
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  INT3(LND) = INT
                  PRINT 11600, LND, IDEF, INT
               END IF
C...channel not implemented
               INT3(8) = 0
               INT3(9) = 0
            END IF
            IF( PARAM .EQ. 'NPVE' ) THEN
               IDEF = NPOVEG
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  NPOVEG = INT
                  PRINT 11700, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'LPIN' ) THEN
C...some options not yet implemented, see manual
               IF( LND .LE. 0 .OR. LND .GT. 11 ) GO TO 30
               IDEF = LPARIN(LND)
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  LPARIN(LND) = INT
                  PRINT 11800, LND, IDEF, INT
               END IF
C...REDEFINITION FOR USE IN NCVS
               LPARIN(12) = 1
               LPAR(1) = 1
               LPAR(2) = LPARIN(2)
               LPAR(3) = LPARIN(3)
               LPAR(4) = LPARIN(1)
               LPAR(7) = LPARIN(7)
               LPAR(8) = LPARIN(8)
               LPAR(9) = LPARIN(9)
               LPAR(10) = LPARIN(10)
               LPAR(11) = 0
               IF( LPARIN(4) .EQ. 1 .OR. LPARIN(5) .EQ. 1 .OR.
     +         LPARIN(6) .EQ. 1 ) LPAR(11) = 1
               LPAR(12) = LPARIN(4)
               LPAR(13) = LPARIN(5)
               LPAR(14) = LPARIN(6)
               LPAR(15) = LPARIN(11)
               LPAR(16) = LPARIN(11)
               LPAR(17) = LPARIN(12)
            END IF
            IF( PARAM .EQ. 'NFMI' ) THEN
               IDEF = NPYMIN
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  IF( INT .GT. 6 ) INT = 6
                  NPYMIN = INT
                  PRINT 11900, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'NFMA' ) THEN
               IDEF = NPYMAX
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  IF( INT .LE. 0 ) INT = 6
                  NPYMAX = INT
                  PRINT 12000, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'ISM2' ) THEN
               IF( LND .LE. 0 .OR. LND .GT. 2 ) GO TO 30
               IDEF = ISAM2(LND)
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  ISAM2(LND) = INT
                  PRINT 12100, LND, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'ISM3' ) THEN
               IF( LND .LE. 0 .OR. LND .GT. 9 ) GO TO 30
               IDEF = ISAM3(LND)
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  ISAM3(LND) = INT
                  PRINT 12200, LND, IDEF, INT
               END IF
               ISAM3(8) = 0
               ISAM3(9) = 0
            END IF
            IF( PARAM .EQ. 'INTO' ) THEN
               IDEF = IOPLOT
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  IOPLOT = INT
                  PRINT 12300, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'WMAS' ) THEN
               RDEF = REAL(MW)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  MW = DBLE(REA)
                  PRINT 12600, RDEF, REA
                  LPAR(5) = 1
               END IF
            END IF
            IF( PARAM .EQ. 'ZMAS' ) THEN
               RDEF = REAL(MZ)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  MZ = DBLE(REA)
                  PRINT 12700, RDEF, REA
                  LPAR(5) = 1
               END IF
            END IF
            IF( PARAM .EQ. 'HMAS' ) THEN
               RDEF = REAL(MH)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  MH = DBLE(REA)
                  PRINT 12800, RDEF, REA
                  LPAR(5) = 1
               END IF
            END IF
            IF( PARAM .EQ. 'TMAS' ) THEN
               RDEF = REAL(MT)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  MT = DBLE(REA)
                  PRINT 12900, RDEF, REA
                  LPAR(5) = 1
               END IF
            END IF
            IF( PARAM .EQ. 'TCIS' ) THEN
               RDEF = REAL(TCUTQ)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  TCUTQ = DBLE(REA)
                  PRINT 13000, RDEF, REA
               END IF
            END IF
            IF( PARAM .EQ. 'TCFS' ) THEN
               RDEF = REAL(TCUTQS)
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  TCUTQS = DBLE(REA)
                  PRINT 13100, RDEF, REA
               END IF
            END IF
            IF( PARAM .EQ. 'IFLO' ) THEN
               IDEF = IFLOPT
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  IFLOPT = INT
                  PRINT 13200, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'M111' ) THEN
               IDEF = MST111
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  MST111 = INT
                  PRINT 13300, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'M115' ) THEN
               IDEF = MST115
               IF( ABS(IDEF-INT) .GT. EPSMAC ) THEN
                  MST115 = INT
                  PRINT 13400, IDEF, INT
               END IF
            END IF
            IF( PARAM .EQ. 'P111' ) THEN
               RDEF = PAR111
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  PAR111 = REA
                  PRINT 13500, RDEF, REA
               END IF
            END IF
            IF( PARAM .EQ. 'P112' ) THEN
               RDEF = PAR112
               IF( ABS(RDEF-REA) .GT. EPSMAC ) THEN
                  PAR112 = REA
                  PRINT 13600, RDEF, REA
               END IF
            END IF
		
		
30      CONTINUE
      ENDIF
*-----------------------------------------------------------------------

10100 FORMAT(' hercha:',
     +'Parameter LST (',I4,') changed from default ',I12,' TO ',I12)
10200 FORMAT(' hercha:',
     +'Parameter PARL(',I4,') changed from default ',E13.6,' TO ',E13.6)
10300 FORMAT(' hercha:',
     +'Parameter EELE  changed from default ',E13.6,' TO ',E13.6)
10400 FORMAT(' hercha:',
     +'Parameter POLA  changed from default ',E13.6,' TO ',E13.6)
10500 FORMAT(' hercha:',
     +'Parameter LLEP  changed from default ',I12,' TO ',I12)
10600 FORMAT(' hercha:',
     +'Parameter EPRO  changed from default ',E13.6,' TO ',E13.6)
10700 FORMAT(' hercha:',
     +'Parameter ICUT  changed from default ',I12,' TO ',I12)
10800 FORMAT(' hercha:',
     +'PARAMETER XMIN  changed from default ',E13.6,' TO ',E13.6)
10900 FORMAT(' hercha:',
     +'Parameter YMIN  changed from default ',E13.6,' TO ',E13.6)
11000 FORMAT(' hercha:',
     +'Parameter YMAX  changed from default ',E13.6,' TO ',E13.6)
11100 FORMAT(' hercha:',
     +'Parameter Q2MI  changed from default ',E13.6,' TO ',E13.6)
11200 FORMAT(' hercha:',
     +'Parameter WMIN  changed from default ',E13.6,' TO ',E13.6)
11300 FORMAT(' hercha:',
     +'Parameter XMAX  changed from default ',E13.6,' TO ',E13.6)
11400 FORMAT(' hercha:',
     +'Parameter EGMN  changed from default ',E13.6,' TO ',E13.6)
11500 FORMAT(' hercha:',
     +'Parameter INT2(',I4,') changed from default ',I12,' TO ',I12)
11600 FORMAT(' hercha:',
     +'Parameter INT3(',I4,') changed from default ',I12,' TO ',I12)
11700 FORMAT(' hercha:',
     +'Parameter NPVE  changed from default ',I12,' TO ',I12)
11800 FORMAT(' hercha:',
     +'Parameter LPIN(',I4,') changed from default ',I12,' TO ',I12)
11900 FORMAT(' hercha:',
     +'Parameter NFMI  changed from default ',I12,' TO ',I12)
12000 FORMAT(' hercha:',
     +'Parameter NFMA  changed from default ',I12,' TO ',I12)
12100 FORMAT(' hercha:',
     +'Parameter ISM2(',I4,') changed from default ',I12,' TO ',I12)
12200 FORMAT(' hercha:',
     +'Parameter ISM3(',I4,') changed from default ',I12,' TO ',I12)
12300 FORMAT(' hercha:',
     +'Parameter INTO  changed from default ',I12,' TO ',I12)
12400 FORMAT(' hercha:',
     +'Parameter MSTA(',I4,') changed from default ',I12,' TO ',I12)
12500 FORMAT(' hercha:',
     +'Parameter PARA(',I4,') changed from default ',E13.6,' TO ',E13.6)
12600 FORMAT(' hercha:',
     +'Parameter WMAS  changed from default ',E13.6,' TO ',E13.6)
12700 FORMAT(' hercha:',
     +'Parameter ZMAS  changed from default ',E13.6,' TO ',E13.6)
12800 FORMAT(' hercha:',
     +'Parameter HMAS  changed from default ',E13.6,' TO ',E13.6)
12900 FORMAT(' hercha:',
     +'Parameter TMAS  changed from default ',E13.6,' TO ',E13.6)
13000 FORMAT(' hercha:',
     +'Parameter TCIS  changed from default ',E13.6,' TO ',E13.6)
13100 FORMAT(' hercha:',
     +'Parameter TCFS  changed from default ',E13.6,' TO ',E13.6)
13200 FORMAT(' hercha:',
     +'Parameter IFLO  changed from default ',I12,' TO ',I12)
13300 FORMAT(' hercha:',
     +'Parameter M111  changed from default ',I12,' TO ',I12)
13400 FORMAT(' hercha:',
     +'Parameter M115  changed from default ',I12,' TO ',I12)
13500 FORMAT(' hercha:',
     +'Parameter P111  changed from default ',E13.6,' TO ',E13.6)
13600 FORMAT(' hercha:',
     +'Parameter P112  changed from default ',E13.6,' TO ',E13.6)
13700 FORMAT(' hercha:',
     +'Parameter THMI  changed from default ',E13.6,' TO ',E13.6)
13800 FORMAT(' hercha:',
     +'Parameter Q2NOT changed from default ',E13.6,' TO ',E13.6)
		
      RETURN
      END
