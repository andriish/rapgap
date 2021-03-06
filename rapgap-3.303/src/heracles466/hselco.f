C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF TS TERM (ELASTIC RADIATIVE TAIL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSELCO(X)
C
C  X(1) -->  XX
C  X(2) -->  G=-1/Q**2   -->  Q**2=-1/G-->  Y
C  X(3) -->  LOG(-TS)
C  X(4) -->  2*K.P
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSKPXY/ XX,Y
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIRCX/ XIRDEL
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSPSPC/ IPHSPC
      DIMENSION X(4)
      COMPLEX*16 HSSRGG,CG
C
C---QUASI-ELASTIC SCATTERING
      XS=1D0
C---CHOOSE VALUE OF Y
      GS=SP-MEI2-MPRO2
      YMAXX=(1D0-4D0*MEI2*MPRO2/GS/GS)/(1D0+2D0*MEI*MPRO/GS)
      IF (ICUT.LT.3) THEN
        GMIN=-1D0/(Q2MIN/XMAX/GS)
        GMAX=-1D0/DMIN1(1D0,YMAXX)
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        GMIN=-1D0/DMAX1(Q2MIN/XMAX/GS,YMIN)
        GMAX=-1D0/DMIN1(1D0,YMAXX,YMAX)
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSELK1'
        STOP
      ENDIF
      GACT=GMIN+X(1)*(GMAX-GMIN)
      Y=-1D0/GACT
      XA=1D0
      CALL HSDELX(XA,Y)
C---X-VALUE
      XXMAX1=HSXMAX(Y)
      XXHH1=1D0-Y-4D0*MEI2*MPRO2/GS/GS
      XXMNY1=(XXHH1+DSQRT(XXHH1*XXHH1-4D0*Y*Y*MEI2*MPRO2/GS/GS))
     &      /2D0/Y/MPRO2*GS
      XXMINY=MEI2/MPRO2/XXMNY1
      XXMIN=DMAX1(XMIN,Q2MIN/Y/GS,XXMINY)
      XXMAX=DMIN1(XMAX,XXMAX1)
      XX=XXMIN+(XXMAX-XXMIN)*X(2)
      Q2=XX*Y*GS
      CALL HSFIVC(XX,Y)
C
      IF(IPRINT.GT.30) THEN
        WRITE(LUNTES,230)SP,XX,Y,XSMIN,XSCUT,DELEPS,DELTA
230     FORMAT(' ***************************************************',/
     F        ,' SP = ',1PD12.3,/
     F        ,' X = ',D12.6,'   Y = ',D12.6,/
     F        ,' XSMIN = ',D17.11,'   XSCUT = ',D17.11,/
     F        ,' DELEPS = ',D12.6,'   DELTA = ',D14.8,/
     F       ,' ***************************************************',//)
      ENDIF
C
      XS=1D0
      CALL HSFCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSELCO=0D0
        RETURN
      ENDIF
      CALL HSFLAB(XX,Y,XS)
      CALL HSLZTS(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1) THEN
        HSELCO=0D0
        RETURN
      ENDIF
      TSMIN=TP+2D0*OMEGA*(ES-EEL-PQ*ZMAX)
      IF(ZMIN.LT.-0.9999D0)THEN
        TSMAX=-XX*XX*XS*MPRO2/(XS-XX+XS*XS*MPRO2/Y/GS)
      ELSE
        TSMAX=TP+2D0*OMEGA*(ES-EEL-PQ*ZMIN)
      ENDIF
C---SUBSTITUTION R = LN(-TS)
      RMAX=DLOG(-TSMIN)
      RMIN=DLOG(-TSMAX)
      R=RMIN+(RMAX-RMIN)*X(3)
      TS=-DEXP(R)
C---NO SUBSTITUTION FOR A1
      CALL HSL1TS(TS,XX,Y,XS,A1MIN,A1MAX,A1M,A1P,CFKP)
      IF(IPHSPC.EQ.1)THEN
       HSELCO=0D0
       RETURN
      ENDIF
      A1=A1MIN+(A1MAX-A1MIN)*X(4)
      SQGRM2=-CFKP*(A1-A1M)*(A1-A1P)
      SQGRAM=DSQRT(DABS(SQGRM2))
      CALL HSFIV1(XX,Y,XS,A1,TS)
C
      A2 = 2D0*DKPS
      RUNALP=1D0
      IF (LPAR(3).GE.3) THEN
        CG=HSSRGG(TS)
        RUNALP=1D0/(1D0+DREAL(CG)/TS)
      ENDIF
      R1=4D0*GX*(
     &     +(T+6D0*MEI2)/(A1-TS)/(A1-TS)/TS
     &     -(T+4D0*MEI2)/(A1-TS)/TS/TS
     &     -(T+2D0*MEI2)/(A2-TS)/(A2-TS)/TS
     &     -(T+4D0*MEI2)/(A2-TS)/TS/TS
     &     +2D0/TS/TS
     &     -1D0/(A1-TS)/TS
     &     +1D0/(A2-TS)/TS
     &     -8D0*MEI2*MEI2/(A1*A2+TS*TS)/TS/TS
     &     +4D0*MEI2*MEI2/(A1*A1+TS*TS)/TS/TS
     &     +4D0*MEI2*MEI2/(A2*A2+TS*TS)/TS/TS )
      FAC1 =-(GS*GS + GU*GU - GX*(GS+GU) -4D0*MEI2*MPRO2)
      FAC2 = -GS*GX + MPRO2*T
      FAC3 = GU*GX - MPRO2*T
      FAC4 = (-2D0*GS*GU + GX*(GS+GU-GX))*2D0*MEI2
      FAC5 = GU*(GX-GU)*2D0*MEI2
      FAC6 = GS*(GX-GS)*2D0*MEI2
      R2=4D0*(
     &      FAC1/(A1+A2-TS)*(1D0/(A1-TS)+1D0/(A2-TS))/TS
     &     -FAC2/(A1-TS)/(A1-TS)/TS
     &     +FAC2/(A1-TS)/TS/TS
     &     -FAC3/(A2-TS)/(A2-TS)/TS
     &     +FAC3/(A2-TS)/TS/TS
     &     +FAC4/(A1*A2+TS*TS)/TS/TS
     &     -2D0*MEI2*MPRO2/(A1-TS)/(A1-TS)/TS
     &     -2D0*MEI2*MPRO2/(A2-TS)/(A2-TS)/TS
     &     -2D0*MPRO2/TS/TS
     &     +MPRO2/(A1-TS)/TS
     &     -MPRO2/(A2-TS)/TS
     &     +FAC5/(A1*A1+TS*TS)/TS/TS
     &     +FAC6/(A2*A2+TS*TS)/TS/TS)
      DO 20 IFL=1,12
 20   CQP(IFL)=0D0
      CALL HSFIE0(-TS,F1EL,F2EL)
      CQP(12)=(F1EL*R1+F2EL*R2)*RUNALP*RUNALP
      SUMME=CQP(12)
      HSELCO=SUMME*Y*2D0*SX1NRM/SQGRAM
     *      *(A1MAX-A1MIN)*(RMAX-RMIN)*(-TS)
     *      *(XXMAX-XXMIN)*(GMAX-GMIN)/GACT**2
      RETURN
      END
