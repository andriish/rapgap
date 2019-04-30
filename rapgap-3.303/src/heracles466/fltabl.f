
C
C***********************************************************************
C
C **********************************************************************
C

      SUBROUTINE FLTABL

C...Integrates the longitudinal structure function, store on grid
C...in x, Q**2.
C   (changed: 12.01.94, HS)
C   (changed: 29.03.95, HS)
chs..use the true value to determine limits of grid also in the case of
c....initial state radiation with reduced effective S stored on PARL(21)

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
C.HS FOR USE IN HERACLES+DJANGO since VERSION 4.5:
C... TRANSFER KINEMATIC LIMITS FROM /HSCUTS/
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      DOUBLE PRECISION MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      DOUBLE PRECISION SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMINH,XMAXH,Q2MINH,Q2MAXH,YMINH,YMAXH,WMINH,GMIN
      DOUBLE PRECISION XMINH,XMAXH,Q2MINH,Q2MAXH,YMINH,YMAXH,WMINH,GMIN
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
C.HS  COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
C.HS &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LINTEG/ NTOT,NPASS
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      EXTERNAL FLQINT,FLGINT,FLTINT
C.HS(12/1/94)
C-HS(10.01.94)
C...Use always fixed center-of-mass energy as upper limit
      PARL21=PARL(21)
      PARL(21)=SNGL(SP-MEI2-MPRO2)
      DO 1 IO=1,15
        IF (INT3(IO).NE.0.OR.ISAM3(IO).NE.0) GOTO 2
    1 CONTINUE
      XMAX=SNGL(XMAXH)
      Q2MIN=SNGL(Q2MINH)
      GOTO 3
    2 CONTINUE
      XMAX=0.999
      Q2MIN=1.0
    3 CONTINUE
      XMIN=SNGL(XMINH)

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      IF(LST(3).GE.3) WRITE(6,1000) LST(11),LQCD,LTM,LHT
      IF(LQCD.LT.1.AND.LTM.LT.1) GOTO 900
      CALL LTIMEX(T1)
      DO 10 IX=1,NFX
      DO 10 IQ=1,NFQ
      FLQT(IX,IQ)=0.
      FLGT(IX,IQ)=0.
   10 FLMT(IX,IQ)=0.
      QR(1)=Q2MIN
      XR(1)=XMIN
      XR(2)=XMAX
      DO 500 IX=1,NFX
      X=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*(IX-1)/(NFX-1))
      QR(2)=X*PARL(21)
      IF(QR(1).GT.QR(2)) GOTO 500
      LQ=0
      DO 400 IQ=1,NFQ
      Q2=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))*
     &(IQ-1)/(NFQ-1))
CTEST IF(LQ.GT.0) GOTO 500
      IF(Q2.GT.PARL(21)) LQ=LQ+1
      Y=Q2/(PARL(21)*X)
      IF(Y.LT.0.0.OR.Y.GT.1.0) LQ=LQ+1
      PARL(25)=ULALPS(Q2)
      IF(LQCD.EQ.1) THEN
C...Quark part.
        ACCUR=PARL(11)
        IT=0
  100   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLQINT,EPS,FLQ)
        IF(FLQ.LT.1) THEN
          ACCUR=FLQ*PARL(11)
          IF(IT.LT.2) GOTO 100
        ENDIF
        FLQT(IX,IQ)=FLQ
C...Gluon part.
        ACCUR=PARL(11)
        IT=0
  200   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLGINT,EPS,FLG)
        IF(FLG.LT.1.) THEN
          ACCUR=FLG*PARL(11)
          IF(IT.LT.2) GOTO 200
        ENDIF
        FLGT(IX,IQ)=FLG
      ENDIF
      IF(LTM.EQ.1) THEN
C...Target mass  part.
        ACCUR=PARL(11)
        IT=0
  300   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLTINT,EPS,FLM)
        IF(FLM.LT.1) THEN
          ACCUR=FLM*PARL(11)
          IF(IT.LT.2) GOTO 300
        ENDIF
        FLMT(IX,IQ)=FLM
      ENDIF
  400 CONTINUE
  500 CONTINUE
  600 CONTINUE
      CALL LTIMEX(T2)
      IF(LST(3).GE.3) WRITE(6,1100) T2-T1
  900 CONTINUE
      PARL(21)=PARL21
      RETURN

 1000 FORMAT(' Initialisation for FL; QCD, target mass, higher twist: ',
     &/,' LST(11) =',I5,' --> LQCD, LTM, LHT =',3I3)
 1100 FORMAT(' FL integrations performed if LQCD=1 and/or LTM=1, ',
     &'results on grid.'/,' Time for FL integrations is ',F7.1,' sec.')
      END
