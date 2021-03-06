C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSSTRF(X,Q2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSNUCL/ HNA,HNZ
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      REAL*4          PYSTOP,PYSLAM
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD
C...HSRADF USED IN SUBROUTINE STRFBS (CALL TO STEIN OR BRASSE)
      COMMON /HSRADC/ AMP , AMP2,RPI ,RPI2 , ALFA, AML , AML2
      COMMON /HSRADF/ AP,AP2,AP4,AL2,AL4,ALPS,CALPI,W2PIT,W2TR,TTR
      DIMENSION DSBOS(2)
      LOGICAL LFIRST
      DATA LFIRST/.TRUE./

      IF (LFIRST) THEN
        LFIRST=.FALSE.
        W2TR=4D0
        TTR=6D0
        AMP2=MPRO2
        AP=2D0*MPRO
        W2PIT=1.15184D0
        W2MIN=1.2321D0
        S=SP
        IF (ILQMOD.EQ.4.AND.ICODE.EQ.0) THEN
          ILIB=2
          ICODE=3034
        ENDIF
        IF (ILQMOD.EQ.5.AND.ICODE.EQ.0) THEN
          ILIB=2
          ICODE=3025
        ENDIF
      ENDIF
      DSBOS(1)=1D0
      DSBOS(2)=Q2/(Q2+MZ2)
      ZF1=0D0
      ZF2=0D0

      IF (ILQMOD.EQ.0.OR.ILQMOD.EQ.1) THEN
C---CALCULATE STRUCTURE FUNCTIONS FROM PARTON DISTRIBUTIONS
        GOTO 1000
      ELSEIF (ILQMOD.EQ.2) THEN
C---FOR LOW Q2 ONLY GAMMA EXCHANGE, PARAMETRIZATION BY BRASSE AND
C   STEIN, COMBINED WITH PARTON DISTRIBUTION FUNCTIONS FOR LARGE Q2
        IF (Q2.LT.TTR) THEN
          DO 1 IB1=1,2
          DO 1 IB2=1,2
          F1(IB1,IB2)=0D0
          F2(IB1,IB2)=0D0
    1     F3(IB1,IB2)=0D0
          CALL STRFBS(X,Q2,ZF1,ZF2)
          F1(1,1)=ZF1
          F2(1,1)=ZF2
          RETURN
         ELSE
          GOTO 1000
        ENDIF
      ELSEIF (ILQMOD.EQ.3.AND.ICODE.NE.0) THEN
C---FOR LOW Q2 ONLY GAMMA EXCHANGE, PARAMETRIZATION BY ALLM (1997) 
C   COMBINED WITH PARTON DISTRIBUTION FUNCTIONS FOR LARGE Q2
C   AND BRASSE PARAMETRIZATION FOR LOW W 
        DO 3 IB1=1,2
        DO 3 IB2=1,2
        F1(IB1,IB2)=0D0
        F2(IB1,IB2)=0D0
    3   F3(IB1,IB2)=0D0
        W2=(1D0-X)/X*Q2+MPRO2
        IF(Q2.LT.TTR) THEN
          IF (W2.LT.W2TR) THEN
            CALL STRFBS(X,Q2,ZF1,ZF2)
            ELSE
            CALL HSSTAL(X,Q2,ZF1,ZF2)
          ENDIF
         ELSEIF (Q2.GE.5000D0.OR.X.GE.0.85D0) THEN
          GOTO 1000
         ELSE
          CALL HSSTAL(X,Q2,ZF1,ZF2)
        ENDIF
        F1(1,1)=ZF1
        F2(1,1)=ZF2
        RETURN
      ELSEIF (ILQMOD.EQ.3.AND.ICODE.EQ.0) THEN
        WRITE(LUNOUT,'(//3(A/))')
     *' ***** WRONG STRUCTURE FUNCTION CODE:'
     *,'       ILQMOD = 3 AND ICODE = 0'
     *,'       EXECUTION STOPPED '
        STOP
      ELSEIF (ILQMOD.EQ.4) THEN
C---BADELEK KWIECINSKI
        DO 4 IB1=1,2
        DO 4 IB2=1,2
        F1(IB1,IB2)=0D0
        F2(IB1,IB2)=0D0
    4   F3(IB1,IB2)=0D0
        IF ((X.GT.0.1D0.AND.Q2.GT.6D0).OR.Q2.GT.998.8D0) THEN
          GOTO 1000
          ELSE
          ANU=Q2/X/2D0
          IF (ANU.LT.10D0) THEN
            CALL STRFBS(X,Q2,ZF1,ZF2)
            ELSE
            CALL HSSTBK(X,Q2,ZF1,ZF2)
          ENDIF
        ENDIF
        F1(1,1)=ZF1
        F2(1,1)=ZF2
        GOTO 2000
      ELSEIF (ILQMOD.EQ.5) THEN
C---DONNACHIE LANDSHOFF with pdf's
        DO 5 IB1=1,2
        DO 5 IB2=1,2
        F1(IB1,IB2)=0D0
        F2(IB1,IB2)=0D0
    5   F3(IB1,IB2)=0D0
        IF (Q2.LT.10D0) THEN
          CALL HSSTDL(X,Q2,ZF1,ZF2)
          F1(1,1)=ZF1
          F2(1,1)=ZF2
          GOTO 2000
         ELSE
          GOTO 1000
        ENDIF
      ELSEIF (ILQMOD.EQ.10) THEN
C---STRUCTURE FUNCTIONS FROM USER ROUTINE
        DO 6 IB1=1,2
        DO 6 IB2=1,2
        F1(IB1,IB2)=0D0
        F2(IB1,IB2)=0D0
 6      F3(IB1,IB2)=0D0
        CALL FIUSER(X,Q2,ZF1,ZF2,IPDFR)
        IF (IPDFR.EQ.1) GOTO 1000
        F1(1,1)=ZF1
        F2(1,1)=ZF2
        GOTO 2000
      ELSE
        WRITE(LUNOUT,200) ILQMOD
        STOP
  200   FORMAT(/,'          WRONG VALUE FOR ILQMOD: ',I8,/
     F          ,' ******** EXECUTION STOPPED IN HSSTRF ',/)
      ENDIF

1000  CONTINUE
      CALL HSPVER(X,Q2)
      IF (HNA.EQ.1D0.AND.HNZ.EQ.1D0) THEN
C---PROTON TARGET
        QUT=QU+QC+QT
        QDT=QD+QS+QB
        QUBT=QBU+QBC+QBT
        QDBT=QBD+QBS+QBB
      ELSEIF (HNA.EQ.1D0.AND.HNZ.EQ.0D0) THEN
C---NEUTRON TARGET
        QDT=QU+QC+QT
        QUT=QD+QS+QB
        QDBT=QBU+QBC+QBT
        QUBT=QBD+QBS+QBB
      ELSE
C---OTHER TARGETS: FIRST CALCULATE DEUTERON STRUCTURE FUNCTIONS
C   PER NUCLEON
        QDT=(QU+QC+QT+QD+QS+QB)/2D0
        QUT=QDT
        QDBT=(QBU+QBC+QBT+QBD+QBS+QBB)/2D0
        QUBT=QDBT
      ENDIF

      DO 10 IB1=1,2
       DO 10 IB2=1,2
         F1(IB1,IB2) = (AFIJ(2,IB1,IB2)*(QUT+QUBT)
     &                 +AFIJ(3,IB1,IB2)*(QDT+QDBT))
     &                 /4D0*DSBOS(IB1)*DSBOS(IB2)  /2D0
         F2(IB1,IB2) = 2D0*X*F1(IB1,IB2)
         F3(IB1,IB2) = (BFIJ(2,IB1,IB2)*(QUT-QUBT)
     &                  +BFIJ(3,IB1,IB2)*(QDT-QDBT))
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)
 10   CONTINUE

C---INCLUDE LONGITUDINAL STRUCTURE FUNCTION
 2000 CONTINUE
      IF (IFLOPT.GT.0) THEN
      FL=0D0
      F2EM=F2(1,1)
      Q2L=Q2
      IF (Q2L.LT.TTR) THEN
        Q2L=TTR
        XP21=SNGL(X)*SNGL(SP-MPRO2-MEI2)
        IF (Q2L.GE.XP21) Q2L=XP21
      ENDIF
      CALL HSLUFL(X,Q2L,F2EM,FL)
      F1(1,1)=(F2(1,1)-FL)/2D0/X
      IF (F1(1,1).LT.0D0) F1(1,1)=0D0
      ENDIF

C---NUCLEAR SHADOWING FOR HEAVY NUCLEI
      IF (HNA.EQ.1D0.AND.HNZ.EQ.1D0) THEN
        CONTINUE
      ELSEIF (HNA.EQ.1D0.AND.HNZ.EQ.0D0) THEN
C--->   CORRECT FOR RATIO F2(NEUTRON)/F2(PROTON)
      ELSE
        DO 11 IB1=1,2
         DO 11 IB2=1,2
         HNRAT=HSNRAT(X)
         F1(IB1,IB2)=F1(IB1,IB2)*HNRAT
         F2(IB1,IB2)=F2(IB1,IB2)*HNRAT
   11   CONTINUE
      ENDIF

C---APPLY LOW Q2 SUPPRESSION
C   MOVED TO LYSTFU
      RETURN
      END
