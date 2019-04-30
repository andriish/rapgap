C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF KP TERM (INITIAL STATE RADIATION)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSTSK1(X)
C
C  X(1) -->  XX
C  X(2) -->  G=-1/Q**2   -->  Q**2=-1/G-->  Y
C  X(3) -->  LOG(XS-XX)
C  X(4) -->  LOG(2*K.P)
C  X(5) -->  TS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSKPXY/ XX,Y
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      COMMON /HSPSPC/ IPHSPC
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSWGTC/ IWEIGS
      DIMENSION X(5),DBOS(2),DSBOS(2),RUNALP(2)
      COMPLEX*16 HSSRGG,CG
C
C---X-VALUE
      XMINL=DLOG(XMIN)
      XMAXL=DLOG(XMAX)
      XXL=XMINL+(XMAXL-XMINL)*X(1)
      XX=DEXP(XXL)
C---Y-VALUE
      GS=SP-MEI2-MPRO2
      YMAXX=XX*(1D0-4D0*MEI2*MPRO2/GS/GS)/(XX*(1D0+XX*MPRO2/GS)+MEI2/GS)
      IF(ICUT.LT.3) THEN
C---CUT IN EXTERNALLY DEFINED Q**2(MIN)
        GMAX=-1D0/(XX*YMAXX*GS)
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        QQ2MIN=XX*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+XX*CTHCON)
        QT2MIN=XX*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=XX*GS/2D0*(1D0-DSQRT(1D0-PTXM0/XX))
        QP2MAX=XX*GS/2D0*(1D0+DSQRT(1D0-PTXM0/XX))
        YP2MAX=QP2MAX/GS/XX
        GMIN=-1D0/MAX(Q2MIN,QQ2MIN,QT2MIN,QP2MIN)
        GMAX=-1D0/(XX*MIN(YMAX,YMAXX,YP2MAX)*GS)
        GMAX=MIN(GMAX,-1D0/Q2MAX)
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSTSK1'
        STOP
      ENDIF
      DG2=DMAX1(GMAX-GMIN,0D0)
C---CUT IN W LATER
      GACT=GMIN+X(2)*DG2
      Y=-1D0/(XX*GACT*GS)
C
C---EXTERNAL WEIGHT
      IACPT=1
      IF (IWEIGS.GT.0) THEN
        CALL HSWGTX(XX,Y,IACPT)
        IF (IACPT.EQ.0) THEN
          HSTSK1=0D0
          RETURN
        ENDIF
      ENDIF
C
      CALL HSFIVC(XX,Y)
      CALL HSDELO(XX,Y)
      XSMIN=HSXSMN(XX,Y)
C     XSCUT=HSXSCT(XX,Y)
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
      XSMAX = 1D0
      XSMINI = XSMIN
      IF(XSMAX.LE.XSMINI) THEN
        HSTSK1=0D0
        RETURN
      ENDIF
C--- SUBSTITUTION UV=LN(XS-X)
      UVMIN = DLOG(XSMINI-XX)
      UVMAX = DLOG(XSMAX-XX)
      UV = UVMIN + (UVMAX-UVMIN)*X(3)
      XS = XX+DEXP(UV)
      CALL HSFCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSTSK1=0D0
        RETURN
      ENDIF
      CALL HSFLAB(XX,Y,XS)
      CALL HSLZK1(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1) THEN
        HSTSK1=0D0
        RETURN
      ENDIF
C---SUBSTITUTION V = LN(A1)
      A1MAX = 2D0*OMEGA*(EEL-PEL*ZMIN)
      IF ((ZMAX.GE.0.9999D0).AND.(EEL/MEI.GT.1D3)) THEN
        A1MIN = 2D0*OMEGA*MEI2/2D0/EEL
      ELSE
        A1MIN = 2D0*OMEGA*(EEL-PEL*ZMAX)
      ENDIF
      VMAX = DLOG(A1MAX)
      VMIN = DLOG(A1MIN)
      V= VMIN + (VMAX-VMIN)*X(4)
      A1 = DEXP(V)
C---NO SUBSTITUTION FOR TS
      CALL HSLTS1(A1,XX,Y,XS,TSMIN,TSMAX,TSM,TSP,CFKP)
      IF(IPHSPC.EQ.1)THEN
       HSTSK1=0D0
       RETURN
      ENDIF
      TS=TSMIN + (TSMAX-TSMIN)*X(5)
      SQGRM2=-CFKP*(TS-TSM)*(TS-TSP)
      IF (SQGRM2.LE.0D0) THEN
        HSTSK1=0D0
        RETURN
      ENDIF
      SQGRAM = DSQRT(DABS(SQGRM2))
      CALL HSFIV1(XX,Y,XS,A1,TS)
C
      A2=2D0*DKPS
      DBOS(1)=1D0
      DBOS(2)=T/(T-MZ2)
      DSBOS(1)=1D0
      RUNALP(1)=1D0
      IF (LPAR(3).GE.3) THEN
        CG=HSSRGG(TS)
        DSBOS(1)=1D0/(1D0+DREAL(CG)/TS)
        RUNALP(1)=DSBOS(1)
      ENDIF
      DSBOS(2)=TS/(TS-MZ2)
      RUNALP(2)=1D0
      IF(LPAR(14).EQ.0)THEN
        SU1=0D0
        SU2=0D0
        PROPI1=0D0
      ELSE
        SU1 = S*S+SS*SS+U*U+US*US
        SU2 = (S-U)*(S+U) + (SS-US)*(SS+US)
        PROPI1= - LLEPT/8D0/XS
     &       * (SS/DKPS/DKQS + S/DKP/DKQ + U/DKPS/DKQ + US/DKP/DKQS)
      ENDIF

      R1=4D0*GX*(
     &     -(T+6D0*MEI2)/(A1-TS)/(A1-TS)/A1
     &     +1D0/(A1-TS)/A1
     &     + 2D0/(A1+A2)/A1
     &     -8D0*MEI2*MEI2/(A1*A2+TS*TS)/(A1+A2)/A1
     &     +4D0*MEI2*MEI2/(A1*A1+TS*TS)/A1/A1
     &     -2D0*MEI2/(A1-TS)/A1/A1 )
C
      FAC1 = -(GS*GS + GU*GU - GX*(GS+GU) -4D0*MEI2*MPRO2)
      FAC2 = -GS*GX + MPRO2*T
      FAC3 = (-2D0*GS*GU + GX*(GS+GU-GX))*2D0*MEI2
      FAC4 = GU*(GX-GU)*2D0*MEI2
      R2=4D0*(
     &     -FAC1/(A1+A2-TS)*(1D0/(A1+A2)+1D0/(A1-TS))/A1
     &     +FAC2/(A1-TS)/(A1-TS)/A1
     &     +FAC3/(A1*A2+TS*TS)/(A1+A2)/A1
     &     -2D0*MPRO2/(A1+A2)/A1
     &     +2D0*MEI2*MPRO2/(A1-TS)/(A1-TS)/A1
     &     +2D0*MEI2*MPRO2/(A1-TS)/A1/A1
     &     -MPRO2/(A1-TS)/A1
     &     +FAC4/(A1*A1+TS*TS)/A1/A1  )
      FK1=-2D0*MEI2*(GS-GU)
      R3= -DFLOAT(LLEPT)*4D0*(
     &     (GS-GU)/(A1+A2)/A1
     &     +FK1/(A1+A2-TS)*(-1D0/(A1+A2)-1D0/(A1-TS))/A1
     &     -(GX/2D0-GS)/(A1-TS)/A1
     &     +GX*T/2D0/(A1-TS)/(A1-TS)/A1
     &     -(GX-2D0*GU)*MEI2/(A1-TS)/(A1-TS)/A1
     &     -(GX-2D0*GU)*MEI2/(A1-TS)/A1/A1)
      DO 20 IFL=1,12
 20     CQP(IFL)=0D0

      SUMME=0D0
      IF (IPDFOP.EQ.0) THEN
C..................................STRUCTURE FUNCTIONS
         CALL HSSTRF(XS,-TS)
         DO 10 IB1=1,2
          DO 10 IB2=1,2
           CQP(12)=CQP(12)+F1(IB1,IB2)*R1*RUNALP(IB1)*RUNALP(IB2)
     &                    +F2(IB1,IB2)*R2*RUNALP(IB1)*RUNALP(IB2)
     &                    +F3(IB1,IB2)*R3*RUNALP(IB1)*RUNALP(IB2)
 10     CONTINUE
        SUMME=CQP(12)
      ELSEIF (IPDFOP.EQ.1) THEN
C..................................PARTON DENSITIES
        CALL HSPVER(XS,-TS)
        AUP=0D0
        AUM=0D0
        ADP=0D0
        ADM=0D0
        FUPI = 0D0
        FDOI = 0D0
        FUPIB =0D0
        FDOIB =0D0
        IF(LPAR(14).NE.0)THEN
          SU3 = (S*S+U*U)/2D0
          SU4 = (S-U)*(S+U)/2D0
          SU5 = (SS*SS+US*US)/2D0
          SU6 = (SS-US)*(SS+US)/2D0
          PPHA1= -T/DKQ/DKQS
          PPHA2= -4D0*MQF2/DKQS/DKQS
          PPHA3= -4D0*MQI2/DKQ/DKQ
          RALEP =
     &     GX/2D0*(
     &     +(T+4D0*MEI2)*(1D0/A2/TS/TS-1D0/A1/TS/TS)
     &     +2D0/TS/TS - 1D0/A1/TS + 1D0/A2/TS + 2D0/A1/A2
     &     -8D0*MEI2*MEI2/A1/A2/TS/TS
     &     +4D0*MEI2*MEI2/TS/TS*(1D0/A1/A1+1D0/A2/A2)
     &     +2D0*MEI2/TS*(1D0/A1/A1+1D0/A2/A2)  )
     &   +XS*(
     &    +FAC1/TS/A1/A2
     &    +FAC2/A1/TS/TS
     &    +(GU*GX - MPRO2*T)/A2/TS/TS
     &    +FAC3/A1/A2/TS/TS
     &    +FAC4/A1/A1/TS/TS
     &    +GS*(GX-GS)*2D0*MEI2/A2/A2/TS/TS
     &    -2D0*MEI2*MPRO2/TS*(1D0/A1/A1+1D0/A2/A2)
     &    -2D0*MPRO2*(1D0/TS/TS+1D0/A1/A2)
     &    +MPRO2/TS*(1D0/A1-1D0/A2) )
          RBLEP =
     &    -DFLOAT(LLEPT)*(
     &     (GS-GU)/A1/A2
     &     -2D0*MEI2*(GS-GU)/A1/A2/TS
     &     +(GX/2D0-GS)/A1/TS
     &     +(GX/2D0-GU)/A2/TS
     &     +GX*T/2D0/A1/TS/TS
     &     +GX*T/2D0/A2/TS/TS
     &     +(GX-2D0*GU)*MEI2/A1/A1/TS
     &     +(2D0*GS-GX)*MEI2/A2/A2/TS)
           PLEP=0D0
           PHAD=0D0
           PFRCH=0D0
      DO 51 IB1=1,2
       DO 51 IB2=1,2
         PLEP=(AFIJ(2,IB1,IB2)*RALEP+BFIJ(2,IB1,IB2)*RBLEP)
     &                *DSBOS(IB1)*DSBOS(IB2)  + PLEP
         PHAD =((AFIJ(2,IB1,IB2)*SU1 - BFIJ(2,IB1,IB2)*SU2*LLEPT)*PPHA1
     &         +(AFIJ(2,IB1,IB2)*SU3 - BFIJ(2,IB1,IB2)*SU4*LLEPT)*PPHA2
     &         +(AFIJ(2,IB1,IB2)*SU5 - BFIJ(2,IB1,IB2)*SU6*LLEPT)*PPHA3)
     &         *DBOS(IB1)*DBOS(IB2)/T/T*4D0/9D0/XS/8D0 +PHAD
         PFRCH=(AFIJ(2,IB1,IB2)*(R1/8D0+R2/4D0*XS)
     &          +BFIJ(2,IB1,IB2)*R3/4D0)
     &                *DSBOS(IB1)*DSBOS(IB2)  + PFRCH
 51      CONTINUE
          SUTOT=PHAD+PLEP
          FRCH=PFRCH/SUTOT
        ELSE
          FRCH=0D0
        ENDIF
        DO 50 IB1=1,2
         DO 50 IB2=1,2
          IF(LPAR(12).NE.0)THEN
           AUP=(AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)+BFIJ(2,IB1,IB2)*R3)
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUP
           AUM=(AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)-BFIJ(2,IB1,IB2)*R3)
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUM
           ADP=(AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)+BFIJ(3,IB1,IB2)*R3)
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADP
           ADM=(AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)-BFIJ(3,IB1,IB2)*R3)
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADM
          ENDIF
         IF(LPAR(14).NE.0)THEN
          PROPI2 = DBOS(IB1)*DSBOS(IB2)/T/TS
          AB12UP = AFIJ(2,IB1,IB2)*SU1 - BFIJ(2,IB1,IB2)*SU2*LLEPT
          AB12UM = AFIJ(2,IB1,IB2)*SU1 + BFIJ(2,IB1,IB2)*SU2*LLEPT
          AB12DP = AFIJ(3,IB1,IB2)*SU1 - BFIJ(3,IB1,IB2)*SU2*LLEPT
          AB12DM = AFIJ(3,IB1,IB2)*SU1 + BFIJ(3,IB1,IB2)*SU2*LLEPT
          FUPI = AB12UP*PROPI1 * PROPI2                   + FUPI
          FDOI = AB12DP*PROPI1 * PROPI2                   + FDOI
          FUPIB =AB12UM*PROPI1 * PROPI2                   + FUPIB
          FDOIB =AB12DM*PROPI1 * PROPI2                   + FDOIB
         ENDIF
 50     CONTINUE
        CQP(1) =  QU *(AUP + FRCH*FUPI  *  2D0/3D0 )
        CQP(2) =  QBU*(AUM + FRCH*FUPIB *(-2D0/3D0)) + CQP(1)
        CQP(3) =  QD *(ADP + FRCH*FDOI  *(-1D0/3D0)) + CQP(2)
        CQP(4) =  QBD*(ADM + FRCH*FDOIB *  1D0/3D0 ) + CQP(3)
        CQP(5) =  QS *(ADP + FRCH*FDOI  *(-1D0/3D0)) + CQP(4)
        CQP(6) =  QBS*(ADM + FRCH*FDOIB *  1D0/3D0 ) + CQP(5)
        CQP(7) =  QC *(AUP + FRCH*FUPI  *  2D0/3D0 ) + CQP(6)
        CQP(8) =  QBC*(AUM + FRCH*FUPIB *(-2D0/3D0)) + CQP(7)
        CQP(9) =  QB *(ADP + FRCH*FDOI  *(-1D0/3D0)) + CQP(8)
        CQP(10) =  QBB*(ADM + FRCH*FDOIB*  1D0/3D0 ) + CQP(9)
        CQP(11) =  QT *(AUP + FRCH*FUPI *  2D0/3D0 ) + CQP(10)
        CQP(12) =  QBT*(AUM + FRCH*FUPIB*(-2D0/3D0)) + CQP(11)
        SUMME=CQP(12)
      ELSEIF (IPDFOP.GE.2) THEN
C..................................PARTON DISTRIBUTION FUNCTIONS
C                                  INCLUDING F_L
        CALL HSSTRF(XS,-TS)
        AUP=0D0
        AUM=0D0
        ADP=0D0
        ADM=0D0
        DO 52 IB1=1,2
         DO 52 IB2=1,2
           AUP=(AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)+BFIJ(2,IB1,IB2)*R3)
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUP
           AUM=(AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)-BFIJ(2,IB1,IB2)*R3)
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUM
           ADP=(AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)+BFIJ(3,IB1,IB2)*R3)
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADP
           ADM=(AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)-BFIJ(3,IB1,IB2)*R3)
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADM
 52     CONTINUE
        CQP(1) =  QU *AUP
        CQP(2) =  QBU*AUM + CQP(1)
        CQP(3) =  QD *ADP + CQP(2)
        CQP(4) =  QBD*ADM + CQP(3)
        CQP(5) =  QS *ADP + CQP(4)
        CQP(6) =  QBS*ADM + CQP(5)
        CQP(7) =  QC *AUP + CQP(6)
        CQP(8) =  QBC*AUM + CQP(7)
        CQP(9) =  QB *ADP + CQP(8)
        CQP(10) = QBB*ADM + CQP(9)
        CQP(11) = QT *AUP + CQP(10)
        CQP(12) = QBT*AUM + CQP(11)
        SUMME=0D0
        DO 53 IB1=1,2
         DO 53 IB2=1,2
          SUMME=SUMME+F1(IB1,IB2)*R1*RUNALP(IB1)*RUNALP(IB2)
     &                +F2(IB1,IB2)*R2*RUNALP(IB1)*RUNALP(IB2)
     &                +F3(IB1,IB2)*R3*RUNALP(IB1)*RUNALP(IB2)
 53     CONTINUE
        RNORM=SUMME/CQP(12)
        DO 54 IF=1,12
 54     CQP(IF)=CQP(IF)*RNORM
      ENDIF

      HSTSK1 = SUMME*Y*2D0*SX1NRM/SQGRAM
     *      * (UVMAX-UVMIN)*(XS-XX)*(VMAX-VMIN)*A1*(TSMAX-TSMIN)
     *      * (XMAXL-XMINL)* XX * DG2/(XX*GS*GACT**2)
      RETURN
      END
