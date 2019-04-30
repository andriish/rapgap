
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE HSPRLG
C---INITIALIZATION OF KINEMATICS / SETTING OF PARAMETERS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIRCX/ XIRDEL
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSNUME/ SIGTOT,SIGTRR,SIGG(20),SIGGRR(20),NEVENT,NEVE(20)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSISGM/ TCUTQ,TCUTQS
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMPLEX*16 CMW2,CMZ2
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSELEP/ IDIPOL
      COMMON /HSNUCL/ HNA,HNZ
      REAL*4          PYSTOP,PYSLAM
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
C----------------------------------
      WRITE(LUNOUT,901)
 901  FORMAT(////,
     1'**************************************************',
     2'*****************************'/)
      WRITE(LUNOUT,'(2(20X,A/))')
     &    'FINAL DEFINITION OF RUN PARAMETERS',
     &    '           FOR HERACLES           '
      WRITE(LUNOUT,902)
 902  FORMAT(
     1'**************************************************',
     2'*****************************')
C
C---ENERGIES AND MOMENTA IN THE HERA LAB SYSTEM
C---ELECTRON BEAM
      IF(EELE.LE.1D0) EELE=3D1
      IF(POLARI.LT.-1D0.OR.POLARI.GT.1D0) POLARI=0D0
      IF(LLEPT.NE.-1.AND.LLEPT.NE.1) LLEPT=-1
C
C---ELECTRO-WEAK PARAMETERS: OPTIONS FOR VIRTUAL&SOFT CONTRIBUTION
C---INCLUDING DEFINITION OF MASSES
      CALL HSSETP
C---TRANSFER TOP-MASS TO PYSTFU
      PYSTOP=MT
C
C---PROTON BEAM
      IF(EPRO.LT.1D0) EPRO=82D1
      PELE=DSQRT((EELE-MEI)*(EELE+MEI))
      PPRO=DSQRT((EPRO-MPRO)*(EPRO+MPRO))
      SP=2D0*(EELE*EPRO+PELE*PPRO)+MPRO2+MEI2
      GSP=SP-MEI2-MPRO2
C
C---KINEMATICAL CUTS
C   DEFINITION OF THE KINEMATICAL CUTS FOR INTEGRATION
C   ICUT=1 : CUTS IN X AND LOWER LIMIT FOR Q2
C       =2 : CUTS IN X AND LOWER LIMITS FOR Q2 AND W
C       =3 : CUTS IN X, Y, AND LOWER LIMITS FOR Q2 AND W
      ICX01=0
      ICX02=0
      ICX03=0
      ICX04=0
      IF (XMIN.LE.0D0) THEN
        XMIN=MEI2/GSP
        ICX01=1
      ENDIF
      IF (XMAX.GE.1D0) THEN
        XMAX=1D0-MEI2/GSP
        ICX02=1
      ENDIF
      IF (YMIN.LE.0D0) THEN
        YMIN=MEI2/GSP
        ICX03=1
      ENDIF
      Q2MNY=MEI2*YMIN*YMIN/(1D0-YMIN)
      IF (Q2MIN.LE.Q2MNY) THEN
        Q2MIN=Q2MNY
        ICX04=1
      ENDIF

      IF(ICUT.LT.1.OR.ICUT.GT.3) ICUT=3
      ICUTO2=0
      ICUTO3=0
      ICUTO4=0
      IF(ICUT.EQ.1) THEN
        IF(Q2MIN.GE.XMIN*GSP) THEN
          XMIN=Q2MIN/GSP
          ICUTO2=1
        ENDIF
        IF (XMIN.GT.XMAX) ICUTO2=4
      ELSEIF(ICUT.EQ.2) THEN
        IF(WMIN.LT.MPRO) WMIN=MPRO
        Q2MIN1=XMIN*(WMIN*WMIN-MPRO2)/(1D0-XMIN)
        XMAX1=1D0-(WMIN*WMIN-MPRO2)/GSP
        IF (Q2MIN1.GT.Q2MIN) THEN
          Q2MIN=Q2MIN1
          ICUTO2=3
        ENDIF
        IF (XMAX1.LT.XMAX) THEN
          XMAX=XMAX1
          ICUTO2=2
        ENDIF
        IF (XMIN.GT.XMAX) ICUTO2=4
      ELSEIF(ICUT.EQ.3) THEN
        IF(WMIN.LT.MPRO) WMIN=MPRO
        YMIN1=Q2MIN/(XMAX*GSP)
        YMIN2=(WMIN*WMIN-MPRO2)/(1D0-XMIN)/GSP
        IF(YMIN1.GT.YMIN.OR.YMIN2.GT.YMIN) ICUTO2=1
        YMIN=MAX(YMIN,YMIN1,YMIN2)
        IF(YMAX.LT.YMIN) ICUTO2=4
        XMIN1=Q2MIN/(YMAX*GSP)
        IF(XMIN1.GT.XMIN.AND.ICUTO2.LE.1) ICUTO2=1
        XMIN=MAX(XMIN,XMIN1)
        XMAX1=1D0-(WMIN*WMIN-MPRO2)/YMAX/GSP
        IF (XMAX1.LT.XMAX) ICUTO2=2
        XMAX=MIN(XMAX,XMAX1)
        IF(XMAX.LT.XMIN) ICUTO2=4
        XMAX2=Q2MAX/YMIN/GSP
        IF (XMAX2.LT.XMAX) ICUTO2=2
        XMAX=MIN(XMAX,XMAX2)
        IF(XMAX.LT.XMIN) ICUTO2=4
        Q2MIN1=XMIN*YMIN*GSP
        Q2MAX1=XMAX*YMAX*GSP
        IF (Q2MIN1.GT.Q2MIN.OR.Q2MAX1.LT.Q2MAX) ICUTO2=3
        Q2MIN=MAX(Q2MIN,Q2MIN1)
        Q2MAX=MIN(Q2MAX,Q2MAX1)
C...CUT ON ELECTRON SCATTERING ANGLE
        IF (THEMIN.NE.0D0) THEN
          CTHMIN=DCOS(THEMIN)
          CTHCON=SP*(1D0+CTHMIN)/4D0/EELE/EELE/(1D0-CTHMIN)
          YMINT1=1D0/(1D0+XMAX*CTHCON)
         ELSE
          CTHMIN=1D0
          CTHCON=1D15
          YMINT1=0D0
        ENDIF
        Q2MINT=XMIN*YMINT1*GSP
        IF (YMIN.LT.YMINT1) THEN
          YMIN=YMINT1
          ICUTO3=1
        ENDIF
        IF (Q2MIN.LE.Q2MINT) THEN
          Q2MIN=Q2MINT
          ICUTO3=1
        ENDIF
C...CUT ON ELECTRON TRANSVERSE MOMENTUM PTMIN
        PTM2=PTMIN*PTMIN
        PTXM0=4D0*PTM2/SP
        YMINP1=(1D0-DSQRT(1D0-PTXM0/XMAX))/2D0
        YMAXP1=(1D0+DSQRT(1D0-PTXM0/XMAX))/2D0
        Q2MPT1=YMINP1*XMAX*GSP
        Q2MPT2=YMAXP1*XMAX*GSP
        IF (XMIN.LT.PTXM0) THEN
          XMIN=PTXM0
          ICUTO4=1
        ENDIF
        XMIN4=0D0
        IF (YMIN.GT.0.5D0) XMIN4=PTM2/SP/YMIN/(1D0-YMIN)
        IF (YMAX.LT.0.5D0) XMIN4=PTM2/SP/YMAX/(1D0-YMAX)
        IF (XMIN.LT.XMIN4) THEN
          XMIN=XMIN4
          ICUTO4=1
        ENDIF
        XMIN5=0D0
        IF (Q2MIN.GT.2D0*PTM2) XMIN5=Q2MIN*Q2MIN/SP/(Q2MIN-PTM2)
        IF (Q2MAX.GT.Q2MPT2) THEN
          Q2MAX=Q2MPT2
          ICUTO4=1
        ENDIF
        IF (Q2MAX.LT.2D0*PTM2) XMIN5=Q2MAX*Q2MAX/SP/(Q2MAX-PTM2)
        IF (XMIN.LT.XMIN5) THEN
          XMIN=XMIN5
          ICUTO4=1
        ENDIF
        IF (YMIN.LT.YMINP1) THEN
          YMIN=YMINP1
          ICUTO4=1
        ENDIF
        IF (YMAX.GT.YMAXP1) THEN
          YMAX=YMAXP1
          ICUTO4=1
        ENDIF
        IF (Q2MIN.LE.Q2MPT1) THEN
          Q2MIN=Q2MPT1
          ICUTO4=1
        ENDIF
      ENDIF
      GMIN=-1D0/Q2MIN
C
C---LOWER LIMIT IN PHOTON ENERGY (OPTIONAL)
      IF (IOPEGM.GT.0) THEN
        DO 2 I=1,5
          INT2(I)=0
          ISAM2(I)=0
  2     CONTINUE
      ENDIF
C
C---PRINT KINEMATICS / BEAM PROPERTIES
C---ELECTRON BEAM
      WRITE(LUNOUT,'(///A/)')
     *    ' *****  PROPERTIES OF THE ELECTRON BEAM  *****'
      WRITE(LUNOUT,'(10X,A,F8.1,A)')
     *      ' ENERGY OF INCIDENT ELECTRON =  ',EELE,' GEV'
      WRITE(LUNOUT,'(10X,A,F8.1,A)')
     *      ' MOMENTUM OF INCIDENT ELECTRON =',PELE,' GEV/C'
      WRITE(LUNOUT,'(10X,A,F11.8,A)')
     *      ' ELECTRON MASS =',ME,' GEV/C**2'
      WRITE(LUNOUT,'(10X,A,I3)')
     *                  ' CHARGE OF INCIDENT ELECTRON =',LLEPT
      WRITE(LUNOUT,'(10X,A,F8.4)')
     *      ' DEGREE OF BEAM POLARIZATION =', POLARI
C
C---PROTON BEAM
      IF (HNA.EQ.1D0.AND.HNZ.EQ.1D0) THEN
        WRITE(LUNOUT,'(///A/)')
     *      ' *****  PROPERTIES OF THE PROTON BEAM  *****'
        WRITE(LUNOUT,'(10X,A,F8.1,A)')
     *      ' ENERGY OF INCIDENT PROTON =  ',EPRO,' GEV'
        WRITE(LUNOUT,'(10X,A,F8.1,A)')
     *      ' MOMENTUM OF INCIDENT PROTON =',PPRO,' GEV/C'
        WRITE(LUNOUT,'(10X,A,F7.4,A)')
     *      ' PROTON MASS =',MPRO,' GEV/C**2'
        WRITE(LUNOUT,'(//10X,A,1PE12.5,A)')
     *      ' CMS ENERGY SQUARED  S =',SP,' GEV**2'
      ELSE
        WRITE(LUNOUT,'(///A/)')
     *      ' *****  PROPERTIES OF THE TARGET BEAM  *****'
        WRITE(LUNOUT,'(10X,A,F4.0,A,F4.0)')
     *      ' A-NUCLEUS = ',HNA,'   Z-NUCLEUS = ',HNZ
        WRITE(LUNOUT,'(10X,A,F8.1,A)')
     *      ' ENERGY PER INCIDENT NUCLEON =  ',EPRO,' GEV'
        WRITE(LUNOUT,'(10X,A,F8.1,A)')
     *      ' MOMENTUM PER INCIDENT NUCLEON =',PPRO,' GEV/C'
        WRITE(LUNOUT,'(10X,A,F7.4,A)')
     *      ' NUCLEON MASS =',MPRO,' GEV/C**2'
        WRITE(LUNOUT,'(//10X,A,1PE12.5,A)')
     *      ' CMS ENERGY SQUARED  S =',SP,' GEV**2'
      ENDIF
C
C---KINEMATICAL CUTS FOR GENERATED EVENTS
      WRITE(LUNOUT,'(///A/)')
     *     ' *****  KINEMATICAL LIMITS FOR GENERATED EVENTS  *****'
      WRITE(LUNOUT,'(10X,A,1PE12.5,2X,A,1PE12.5)')
     *     ' XMIN=',XMIN,' XMAX=',XMAX
      WRITE(LUNOUT,'(10X,A,1PE12.5,A)')
     *     ' Q2MIN=',Q2MIN,' GEV**2, '
      WRITE(LUNOUT,'(10X,A,1PE12.5,A)')
     *     ' Q2MAX=',Q2MAX,' GEV**2'
      WRITE(LUNOUT,'(10X,A,1PE12.5,A,17X,A)')
     *     ' WMIN=', WMIN, ' GEV',' (ACTIVE ONLY FOR ICUT=2)'
      WRITE(LUNOUT,'(10X,A,1PE12.5,2X,A,1PE12.5,1X,A)')
     *     ' YMIN=',YMIN,' YMAX=',YMAX, ' (ACTIVE ONLY FOR ICUT=3)'
      WRITE(LUNOUT,'(10X,A,I2)')
     *     ' ICUT=',ICUT
C
      IF(ICUTO2.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' NOTE: XMIN AND/OR YMIN MODIFIED FOR CONSISTENCY WITH'
        WRITE(LUNOUT,'(10X,A)')
     &     ' MINIMUM ALLOWED Q**2 OR MINIMUM ALLOWED W'
      ELSEIF(ICUTO2.EQ.2) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' NOTE: XMAX MODIFIED FOR CONSISTENCY WITH'
        WRITE(LUNOUT,'(10X,A)')
     &     ' MIN/MAX ALLOWED Q**2 AND MAX/MIN ALLOWED Y'
      ELSEIF(ICUTO2.EQ.3) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' NOTE: Q2MIN/Q2MAX MODIFIED FOR CONSISTENCY WITH'
        WRITE(LUNOUT,'(10X,A)')
     &     ' MINIMUM ALLOWED W, X, Y'
      ELSEIF(ICUTO2.EQ.4) THEN
        WRITE(LUNOUT,'(/10X,A/10X,A/10X,A)')
     &     ' NOTE: INCONSISTENT KINEMATICAL LIMITS: ',
     &     '       XMIN > XMAX AND/OR YMIN > YMAX ',
     &     '       AFTER CONSISTENCY CHECK QITH Q**2_MIN AND W_MIN. ',
     &     '       EXECUTION STOPPED'
        STOP
      ENDIF
      IF (ICUTO3.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,2A)')
     &     ' NOTE: YMIN AND/OR Q2MIN MODIFIED FOR CONSISTENCY WITH',
     &     ' CUT ON ELECTRON SCATTERING ANGLE'
      ENDIF
      IF (ICUTO4.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,2A)')
     &     ' NOTE: YMIN/MAX, XMIN AND/OR Q2MIN/MAX MODIFIED FOR',
     &     ' CONSISTENCY WITH CUT ON ELECTRON TRANSVERSE MOMENTUM'
      ENDIF
      IF (ICX01.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' XMIN MODIFIED, VALUES <= 0 NOT ALLOWED '
      ENDIF
      IF (ICX02.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' XMAX MODIFIED, VALUES >= 1 NOT ALLOWED '
      ENDIF
      IF (ICX03.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' YMIN MODIFIED, VALUES <= 0 NOT ALLOWED '
      ENDIF
      IF (ICX04.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' Q2MIN MODIFIED, VALUES <= 0 NOT ALLOWED '
      ENDIF
C
C---LOWER LIMIT IN PHOTON ENERGY (OPTIONAL)
      IF (IOPEGM.GT.0) THEN
        WRITE(LUNOUT,'(/10X,A,1PE10.3,A)')
     &       ' MINIMUM PHOTON ENERGY REQUESTED: EGMIN = ',EGMIN,' GEV'
        WRITE(LUNOUT,'(10X,A)')
     &       ' (EVENT SAMPLING FOR HARD-PHOTON CONTRIBUTIONS ONLY)'
      ENDIF
C
C---PARAMETERS FOR GSW THEORY
      WRITE(LUNOUT,'(///A/)')
     *    ' *****  PARAMETERS FOR EL-WEAK THEORY  *****'
      WRITE(LUNOUT,'(10X,A,I2,10X,A/29X,A/)')
     *    ' IFIX =', LPAR(4), ' IFIX=1 : MW FIXED',
     *                        '     =2 : GF FIXED'
      WRITE(LUNOUT,212) MW,MZ,SW2
      WWIDTH=-DIMAG(CMW2)/DSQRT(DREAL(CMW2))
      ZWIDTH=-DIMAG(CMZ2)/DSQRT(DREAL(CMZ2))
      WRITE(LUNOUT,213) WWIDTH,ZWIDTH
      WRITE(LUNOUT,214) MU,MC,MD,MB,MS,MT
      WRITE(LUNOUT,215) MH
212   FORMAT(10X,' BOSON MASSES:   MW = ',F10.4,' GEV,',/
     F       10X,'                 MZ = ',F10.4,' GEV,   SW2 = ',F10.4)
213   FORMAT(10X,' BOSON WIDTHS:   GW = ',F10.4,' GEV,    GZ = ',F10.4)
214   FORMAT(10X,
     F       ' FERMION MASSES: MU = ',F10.4,' GEV,  MC = ',F10.4,' GEV',
     F /,10X,'                 MD = ',F10.4,' GEV,  MB = ',F10.4,' GEV',
     F /,10X,'                 MS = ',F10.4,' GEV,  MT = ',F10.4,' GEV')
215   FORMAT(10X,' HIGGS MASS:     MH = ',F10.4,' GEV')
C
C---PARTON DISTRIBUTION
      WRITE(LUNOUT,'(///A/A/)')
     *    ' *****  OPTIONS FOR PARTON DISTRIBUTIONS OR      *****',
     *    ' *****  STRUCTURE FUNCTIONS                      *****'
      CALL HSWPDF
C
C---ELASTIC SCATTERING
      IF (IDIPOL.NE.0) THEN
      WRITE(LUNOUT,'(///A/A/)')
     *    ' *****  ELASTIC SCATTERING INCLUDED               *****',
     *    ' *****  WITH DIPOLE FORM FACTOR FROM STEIN ET AL. *****'
      ENDIF
C
C---ANGULAR CUTS FOR QUARKONIC BREMSSTRAHLUNG
      IF(INT3(4).GE.1 .OR. ISAM3(4).GE.1 .OR.
     &   ((INT2(1).GE.1.OR.ISAM2(1).GE.1).AND.LPARIN(5).GE.1) ) THEN
        WRITE(LUNOUT,'(//A/)')
     &          ' *****  ANGULAR CUTS FOR QUARKONIC RADIATION  *****'
        WRITE(LUNOUT,'(10X,A,1PE10.3,A)')
     &    ' TCUTQ  =', TCUTQ, ' RAD',
     &    ' TCUTQS =', TCUTQS, ' RAD'
      ENDIF
C
C---OPTIONS FOR VIRTUAL&SOFT CONTRIBUTION
      IF(INT2(1).EQ.1.OR.INT2(2).EQ.1
     *   .OR.ISAM2(1).GE.1.OR.ISAM2(2).GE.1) THEN
        WRITE(LUNOUT,'(//A/)')
     *    ' *****  OPTIONS FOR VIRT&SOFT CONTRIBUTION         *****'
        WRITE(LUNOUT,216) LPAR
216     FORMAT(/,' PARAMETER LIST',/
     F,' *******************************************************',/
     F,' **      BORN CROSS SECTION:    LPAR( 1) = ',I6,'     **',/
     F,' **      1-LOOP CORRECTIONS:    LPAR( 2) = ',I6,'     **',/
     F,' **      HIGHER ORDERS:         LPAR( 3) = ',I6,'     **',/
     F,' **      MW OR GMU FIXED:       LPAR( 4) = ',I6,'     **',/
     F,' **      MASS PARAM FROM INPUT: LPAR( 5) = ',I6,'     **',/
     F,' **      STRUCTURE FUNCTIONS:   LPAR( 6) = ',I6,'     **',/
     F,' **      SIGMA-GAMMA:           LPAR( 7) = ',I6,'     **',/
     F,' **      SIGMA-GAMMA-Z:         LPAR( 8) = ',I6,'     **',/
     F,' **      SIGMA-Z:               LPAR( 9) = ',I6,'     **',/
     F,' **      SIGMA-W:               LPAR(10) = ',I6,'     **',/
     F,' **      QED CORRECTIONS:       LPAR(11) = ',I6,'     **',/
     F,' **      LEPTONIC QED:          LPAR(12) = ',I6,'     **',/
     F,' **      HADRONIC QED:          LPAR(13) = ',I6,'     **',/
     F,' **      LEPTON-QUARK-INTRF:    LPAR(14) = ',I6,'     **',/
     F,' **      WEAK CORRECTIONS:      LPAR(15) = ',I6,'     **',/
     F,' **      WEAK BOXES:            LPAR(16) = ',I6,'     **',/
     F,' **      GAMMA OR/AND Z         LPAR(17) = ',I6,'     **',/
     F,' **      NOT USED:              LPAR(18) = ',I6,'     **',/
     F,' **      NOT USED:              LPAR(19) = ',I6,'     **',/
     F,' **      NOT USED:              LPAR(20) = ',I6,'     **',/
     F,' *******************************************************'/)
      ENDIF
C
C---ELASTIC CONTRIBUTIONS COMPATIBLE WITH KINEMATIC CUTS?
      IF (XMAX.LT.(1D0-MEI2/GSP-XIRDEL)
     &   .AND.(INT2(3).NE.0.OR.ISAM2(3).NE.0))THEN
        INT2(3)=0
        ISAM2(3)=0
        WRITE(LUNOUT,'(/A,3(/A,10X,A))')
     &    ' ********** WARNING **********',
     &    ' ***',' VALUE OF XMAX DOES NOT ALLOW NON-RADIATIVE ',
     &    ' ***',' ELASTIC EP SCATTERING ',
     &    ' ***',' INT2(3) AND ISAM2(3) SET TO 0 '
      ENDIF
      IF (WMIN.GT.MPRO.AND.(INT2(3).NE.0.OR.ISAM2(3).NE.0.OR.
     &   INT3(10).NE.0.OR.INT3(11).NE.0.OR.INT3(12).NE.0.OR.
     &   ISAM3(10).NE.0.OR.ISAM3(11).NE.0.OR.ISAM3(12).NE.0)) THEN
        INT2(3)=0
        INT3(10)=0
        INT3(11)=0
        INT3(12)=0
        ISAM2(3)=0
        ISAM3(10)=0
        ISAM3(11)=0
        ISAM3(12)=0
        WRITE(LUNOUT,'(/A,3(/A,10X,A))')
     &    ' ********** WARNING **********',
     &    ' ***',' VALUE OF WMIN DOES NOT ALLOW ',
     &    ' ***',' ELASTIC AND QUASI-ELASTIC EP SCATTERING ',
     &    ' ***',' INT2(3), ISAM2(3), INT3(10-12) AND ISAM3(10-12) SET
     &TO 0 '
      ENDIF
C
C---OPTIONS FOR INTEGRATION / SAMPLING
      WRITE(LUNOUT,'(///A/)')
     *    ' *****  OPTIONS FOR INTEGRATION / SAMPLING  *****'
C---CHANNELS FOR CHARGED CURRENT NOT YET DEFINED
      IF (INT3(8).NE.0.OR.ISAM3(8).NE.0) THEN
        INT3(8)=0
        ISAM3(8)=0
        LPAR(13)=0
        WRITE(LUNOUT,'(10X,A,/,10X,A)')
     *  ' CHARGED CURRENT: QUARKONIC RADIATION NOT YET ACTIVATED ',
     *  ' ICC32, ISCC32 AND LPARIN(5) SET TO ZERO ' 
      ENDIF
      IF (INT3(9).NE.0.OR.ISAM3(9).NE.0) THEN
        INT3(9)=0
        ISAM3(9)=0
        LPAR(14)=0
        WRITE(LUNOUT,'(10X,A,/,10X,A)')
     *' CHARGED CURRENT: LEPTON-QUARK INTERFERENCE NOT YET ACTIVATED ',
     *' ICC33, ISCC33 AND LPARIN(6) SET TO ZERO ' 
      ENDIF

      DO 1 I=1,5
        IF(INT2(I).LT.0.OR.INT2(I).GT.1) INT2(I)=0
 1    CONTINUE
      WRITE(LUNOUT,'(10X,A,8I4)')
     *      ' INT2(I) ', INT2
      WRITE(LUNOUT,'(10X,A,15I4)')
     *      ' INT3(I) ', INT3
      WRITE(LUNOUT,'(10X,A,8I4)')
     *      ' ISAM2(I)',ISAM2
      WRITE(LUNOUT,'(10X,A,15I4)')
     *      ' ISAM3(I)', ISAM3
      IF(LPAR(14).GT.0) THEN
        ISM3TT=ISAM3(1)+ISAM3(2)+ISAM3(3)+ISAM3(4)
        IF ((ISM3TT.GT.0).AND.(ISAM3(1).EQ.0.OR.ISAM3(2).EQ.0.OR.
     *                         ISAM3(3).EQ.0.OR.ISAM3(4).EQ.0)   ) THEN
        WRITE(LUNOUT,'(/A,3(/A,10X,A))')
     &    ' ********** WARNING **********',
     &    ' ***',' LEPTON-QUARK INTERFERENCE CAN BE INCLUDED ONLY IF',
     &    ' ***',' AT THE SAME TIME ALL LEPTONIC AND QUARKONIC ',
     &    ' ***',' CHANNELS ARE REQUESTED,  LPARIN(6) SET TO 0 '
        LPARIN(6)=0
        LPAR(14)=0
        ENDIF
      ENDIF
C
C---NUMBER OF REQUESTED EVENTS
      IF (ISAM2(1).NE.0 .OR. ISAM2(2).NE.0.OR.ISAM2(3).NE.0
     &    .OR. ISAM3(1).GT.0 .OR. ISAM3(2).GT.0 .OR. ISAM3(3).GT.0
     &    .OR. ISAM3(4).GT.0 .OR. ISAM3(7).GT.0 .OR. ISAM3(9).GT.0
     &    .OR. ISAM3(10).GT.0 .OR. ISAM3(11).GT.0 .OR. ISAM3(12).GT.0
     &   ) THEN
        WRITE(LUNOUT,'(///A/)')
     *              ' *****  NUMBER OF EVENTS TO BE SAMPLED  *****'
        WRITE(LUNOUT,'(10X,A,I8)')
     *    ' NUMBER OF EVENTS REQUESTED  NEVENT =',NEVENT
      ENDIF
C---
      WRITE(LUNOUT,'(///)')
      RETURN
      END
