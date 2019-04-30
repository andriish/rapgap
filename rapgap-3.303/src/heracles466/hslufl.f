C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
************************************************************************
C
C   SUBROUTINES FOR THE INCLUSION OF THE LONGITUDINAL STRUCTURE
C   FUNCTION IN HERACLES, since VERSION 4.5
C
C   TAKEN / MODIFIED FROM LEPTO 6.1 BY G.INGELMAN
C
************************************************************************

      SUBROUTINE HSLUFL(XA,Q2A,F2EM,FL)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C...COMMON BLOCKS FROM HERACLES (NOTE: NAMES HAVE PARTLY BEEN CHANGED)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCUTS/ XMINH,XMAXH,Q2MINH,Q2MAXH,YMINH,YMAXH,WMINH,GMIN
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      REAL*4          PYSTOP,PYSLAM
C...COMMON BLOCKS FROM LEPTO
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      REAL*4          CUT            ,PARL    ,X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      REAL*4          PARI    ,EWQC       ,QC   ,ZL     ,ZQ     ,PQ
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      REAL*4                  XR   ,QR   ,FLQT       ,FLGT       ,
     &FLMT
      REAL*4 FLQ,FLG,FLM
      LOGICAL LFIRST,LFIRS1
      DATA LFIRST,LFIRS1 /2*.TRUE./

      IF (LFIRS1) THEN
      LFIRS1=.FALSE.
      QC(1)=-.33333
      QC(2)=.66667
      QC(3)=-.33333
      QC(4)=.66667
      QC(5)=-.33333
      QC(6)=.66667
      QC(7)=-.33333
      QC(8)=.66667
C...ONLY WARNINGS PRINTED, EXECUTION NOT STOPPED IF ERROR IN LEPTO
ckc..
      LST3=LST(3)
      LST(3)=1
C...F_L INCLUDED
      LST(11)=IFLOPT
C...To be sure that LNSTRF is initialized
      CALL HSPVER(0.1D0,100D0)
      ENDIF

      FLQ=0D0
      FLG=0D0
      FLM=0D0
      FLT=0D0
      X=SNGL(XA)
      Q2=SNGL(Q2A)
      IF (LQCD.EQ.1.OR.LTM.EQ.1) THEN
C...F_L FROM INTERPOLATION
        IF (LFIRST) THEN
C...INITIALIZATION OF GRID
        LFIRST=.FALSE.
        PARL(21)=SNGL(SP-MEI2-MPRO2)
        NFX=41
        NFQ=16
        DO 1 I=1,NFX
        DO 1 J=1,NFQ
        FLQT(I,J)=0.
    1   FLGT(I,J)=0.
        CALL FLTABL
        ENDIF
        CALL FLIPOL(FLQ,FLG,FLM)
      ELSEIF (LQCD.EQ.2.OR.LTM.EQ.2) THEN
C...F_L FROM INTEGRATION EVENT-BY-EVENT
        IF (LFIRST) THEN
C...INITIALIZATION
        LFIRST=.FALSE.
C...PROTON TARGET
        PARL(1)=1.
        PARL(2)=1.
        PARI(11)=(PARL(1)-PARL(2))/PARL(1)
C...KINEMATIC LIMITS
        XR(1)=SNGL(XMINH)
        XR(2)=SNGL(XMAXH)
        QR(1)=SNGL(Q2MINH)
        ENDIF
        CALL FLINTG(FLQ,FLG,FLM)
      ELSE
        RETURN
      ENDIF

      IF (LHT.GE.1) FLT=8D0*PARL(19)/Q2A*F2EM
      FL=DBLE(FLQ+FLG+FLM)+FLT
      LST(3)=LST3
      RETURN
      END
