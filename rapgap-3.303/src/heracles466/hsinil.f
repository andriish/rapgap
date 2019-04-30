C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C***********************************************************************
C
C     NEW SUBROUTINES FOR THE INCLUSION OF
C     ELASTIC AND QUASI-ELASTIC EP SCATTERING
C     AND NUCLEAR SHADOWING
C
C
C     NEW ROUTINES:      IN HSOURCEE: HSELG1 (F)
C                                     HSELG2 (F)
C                                     HSEL22 (F)
C                                     HSELK1 (F)
C                                     HSELK2 (F)
C                                     HSELCO (F)
C                                     HSSGEL (F)
C                                     HSXMAX (F)
C                                     HSNRAT (F)
C                                     HSINIL (S)
C                                     HSFIE0 (S)
C                                     HSFIEL (S)
C                                     HSDELX (S)
C                                     D01AJF (S)
C
C***********************************************************************
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSINIL(FUN,EPSO,NBIN2,NDO2,SIG2,SIG2E,XX2)
C
C   INITIALIZATION FOR ELASTIC EP EVENT GENERATION
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C      INTEGER MINPTS,MAXPTS
      EXTERNAL FUN
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSINTL/ XL,XU
      DIMENSION XX2(50,1)
      DATA BLOW/0D0/, BUP/1D0/
C      DATA MINPTS,MAXPTS / 0, 1000/
      PARAMETER (LW=2000,LIW=500)
      DIMENSION IW(LIW),W(LW)
C
C   PARAMETERS
C
      ACCURA=EPSO
      NDO2=NBIN2
      SIG2=0D0
      SIG2E=0D0
      IFAIL=1
      NDIMEN=1
      CALL D01AJF(FUN,BLOW,BUP,0D0,ACCURA,RESULT,ACCFIN,
     *                W,LW,IW,LIW,IFAIL)
      IF(IFAIL.NE.0.OR.IPRINT.GT.1)
     &  WRITE(LUNTES,'(A)')
     &  ' D01AJF DID NOT MEET REQUIRED ACCURACY IN ELASTIC EP '
      SIG2=RESULT
      SIG2E=ACCFIN
C
      IF(IPRINT.GT.1) THEN
        WRITE(LUNOUT,'(///A,5X,1PE12.4,A,1PE12.4,A)')
     *        ' CROSS SECTION VALUE SIG2L (WITH ERROR ESTIMATE):',
     *        SIG2, ' +/- ', SIG2E, '  NB'
        WRITE(LUNTES,'(A,1PD10.1)') ' RELATIVE ACCURACY REQUIRED:',
     &        ACCURA
      ENDIF
      G=0D0
      DG=1D0/DFLOAT(NDO2)
      DO 5 I=1,NDO2
        XX2(I,1)=DFLOAT(I)*DG
  5   CONTINUE
      IF(IPRINT.GT.1) THEN
        WRITE(LUNTES,'(A,/,4(5(1PD15.5)/))')
     &                     ' XX2(I)', (XX2(IG,1),IG=1,NDO2)
      ENDIF
      IF(IPRINT.GT.1) WRITE(LUNTES,'(A)') ' HSINIL FINISHED'
      RETURN
      END
