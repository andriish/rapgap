
      FUNCTION PQCDQQB(X1,Y1)
      IMPLICIT None
	Double Precision X,XV,WMAX
	Integer Ndimen,IMIX
      DIMENSION X(20)
      REAL X1,Y1,PQCDQQB
      COMMON /XVAL/ XV(20),NDIMEN
      COMMON /OALPHAS/ WMAX,IMIX
#include "rgdisdif.inc"
      REAL SNGL
      Double Precision draprn
	Integer I
	Double Precision WDUM,WMAX1,FXN1
      EXTERNAL draprn
      DATA WDUM/0./
      IF(IDIR.EQ.0) THEN
C generate additional random numbers for order alpha_s QCD processes
C           X(NDIMEN+1) = phi in PHASE
C           X(NDIMEN+2) = cost in PHASE
C           X(NDIMEN+3) = XP2 = E_part/E_proton
C           X(NIDMEN+4) = phi electron
         XV(NDIMEN+1) = draprn()
         XV(NDIMEN+4) = draprn()
         XV(NDIMEN+2) = DBLE(X1)
         XV(NDIMEN+3) = DBLE(Y1)
      ELSEIF(IDIR.EQ.1) THEN
C generate additional random numbers for order alpha_s QCD processes
C           X(NDIMEN+1) = XP2 = E_part/E_proton
C           X(NIDMEN+2) = phi electron
C           X(NDIMEN+3) = phi in PHASE
C           X(NDIMEN+4) = cost in PHASE
         XV(NDIMEN+1) = DBLE(X1)
         XV(NDIMEN+2) = draprn()
         XV(NDIMEN+3) = draprn()
         XV(NDIMEN+4) = DBLE(Y1)
      ELSE
         WRITE(6,*) ' PQCDQQB: selection not possible. IDIR = ',IDIR
      ENDIF
c      write(6,*) (XV(I),I=Ndimen+1,Ndimen+4)
      DO 10 I=1,NDIMEN+4
   10 X(I) = XV(I)
c      write(6,*) ' pqcdqqb: ndimen ',NDIMEN,IDIR
c      write(6,*) ' pqcdqqb: x(i) ',(X(I),I=1,NDIMEN+4)
      if(x1.gt.1..or.x1.lt.0..or.y1.gt.1..or.y1.lt.0.) THEN
         write(6,*) ' x1 = ',x1,' y1 = ',y1
      endif
      WMAX1 = FXN1(X,WDUM)
c      write(6,*) ' pcdqqb: fxn1 wmax',wmax1,wmax,idir,ipro
      IF(WMAX1.GT.WMAX) WMAX=WMAX1
      PQCDQQB = SNGL(WMAX1)
      RETURN
      END
      FUNCTION FL(X)
	REAL FL
      FL=0.
      RETURN
      END
      FUNCTION FU(X)
	REAL FU
      FU=1.
      RETURN
      END
