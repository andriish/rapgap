C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSNC22(X,Q2)
C
C   D2SIG/DX*DQ2 WITH COMPLETE 1-LOOP SOFT AND VIRTUAL CORRECTIONS
C   FROM H. SPIESBERGER
C
C        NOTE: OUTPUT FROM H.S.  DSIG / DX*DY
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
C---------------------------------------------------------------------
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKPXY/ XX,Y
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARM/ POLARI,LLEPT,LQUA
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      DATA NEVERR /0/
      LOGICAL LERR
      DATA LERR /.TRUE./
C
      XX=X
      IF(IPRINT.GT.20)
     &  WRITE(LUNTES,'(A/3(1PD13.5),F8.3,2I3)')
     *         ' HSNC22: SP, X, Q2, POLARI,LLEPT,LQUA',
     *         SP,X,Q2,POLARI,LLEPT,LQUA
      Y=Q2/X/(SP-MEI2-MPRO2)
      HSNC22=HSSGNC(X,Y,LLEPT,POLARI,LQUA)/X/SP
C
      IF(HSNC22.LT.0D0) THEN
        HSNC22=0D0
        NEVERR=NEVERR+1
        IF (NEVERR.LT.10) THEN
         WRITE(LUNTES,'(A,/,4(1PD13.5),2I3,F8.3/A/2(6(1PD13.5)/))')
     +     ' HSNC22: X, Y, Q2, HSNC22, LLEPT, LQUA, POLARI',
     +      X, Y, Q2, HSNC22, LLEPT, LQUA, POLARI,
     +     '     HSPDFQ: QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT',
     +      QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
        ELSEIF (LERR) THEN
         LERR=.FALSE.
         WRITE(LUNTES,'(A,I3,A)')
     &     ' ERROR HSNC22 < 0 HAS OCCURED ',NEVERR,
     &     ' TIMES, NO FURTHER WARNINGS ARE PRINTED'
        ELSE
        ENDIF
      ENDIF
      RETURN
      END
