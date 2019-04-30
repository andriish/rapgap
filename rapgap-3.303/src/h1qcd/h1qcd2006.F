
      SUBROUTINE H1QCD2006(BETA_IN,Q2_IN,XPQ,X_POM_IN,T2_IN)

      IMPLICIT REAL*8 (A-G,O-Z)
      REAL*4 T2_IN,X_POM_IN,BETA_IN,Q2_IN,XPQ(-6:6)
      DIMENSION XPQP(-6:6),XPQM(-6:6)
      double precision F2(2),FL(2),C2(2),CL(2)
      INTEGER ICALL
      DATA ICALL/0/
      save icall
      character com*72

#include "rgdiffr.inc"
#include "rgnlodiff.inc"
#include "rgpdiss.inc"

      icall = icall + 1

C     Set Parameters for fits

c      WRITE(6,*)'H1QCD2006',NPOM,NG

      IF(icall.eq.1) THEN
c         open(unit=70,FILE=".message",STATUS="NEW")
c         write(70,*) 'The precious Diffractive H1 Fit 2002 is used '
c         write(70,*) ' '
c         close(70)
c         com="mail jung@mail.desy.de  -s THE_H1_FIT_IS_USED < .message"
c         call system(com)
c         com="rm .message"
c         call system(com)
         WRITE(6,*)'#############################################'
         WRITE(6,*)'#           H1QCD 2006 fits Selected        #'
         WRITE(6,*)'#############################################'
         WRITE(6,*)'# FIT is ',NG
         if ((abs(ng)-29).eq.1) write(6,*)'# fit A '
         if ((abs(ng)-29).eq.2) write(6,*)'# fit B '
         IF(NPOM.EQ.-30) WRITE(6,*)'# Pomeron Part Only'
         IF(NPOM.EQ.-31) WRITE(6,*)'# Reggeon Part Only'
         IF(NPOM.EQ.-32) WRITE(6,*)'# Pomeron + Reggeon'
         WRITE(6,*)'#############################################'
      ENDIF

C     Input quantities are REAL*4 -> Convert to REAL*8 ------------------------

      Z =DBLE(BETA_IN)
      Q2   =DBLE(Q2_IN)
      X_POM=DBLE(X_POM_IN)
      T2   =DBLE(T2_IN)

c     Reggeon Init pdf's
c     PDFLIB-Owens(2,1,1)
      IF (NPOM.EQ.-31 .OR. NPOM.EQ.-32) THEN
c        write(6,*) ' LHApdf called '
        CALL PION_INIT(0.D0,0.D0,211.D0) ! LO-Owens Set 1 
cc        write(6,*) ' pdflib called '
c        CALL PION_INIT(2.D0,1.D0,1.D0) ! LO-Owens Set 1
      ENDIF

C     Final DIS 2006 Fits from paper
      IF(NPOM.LE.-30 .AND. NPOM.GE.-32) THEN
         IF(icall.eq.1) THEN
            IFIT=ABS(NG)-29
         ELSE
            IFIT=0
         ENDIF
c pomeron pdfs
         CALL QCD_2006(Z,Q2,IFIT,XPQP,F2,FL,C2,CL)
c reggeon pdfs
         IF (NPOM.EQ.-31 .OR. NPOM.EQ.-32) CALL pion_xpq(Z,Q2,XPQM)
c pomeron flux
         IFFIT=ABS(NG)-29
         call h12006flux(x_pom,t2,0,iffit,1,flux_p)
c reggeon flux
         call h12006flux(x_pom,t2,0,iffit,2,flux_m)

      ENDIF

      DO I=-6,6
         XPQ(I)=0.
         IF(ABS(I).LT.6) THEN
            IF(NPOM.EQ.-30) THEN
               XPQ(I)=REAL(XPQP(I)*FLUX_P)
		   Diff_flux = FLUX_P
            ELSEIF(NPOM.EQ.-31 ) THEN
               XPQ(I)=REAL(XPQM(I)*FLUX_M)
		   Diff_flux = FLUX_M		   
            ELSEIF(NPOM.EQ.-32) THEN
               XPQ(I)=REAL(XPQM(I)*FLUX_M+XPQP(I)*FLUX_P)
		   Diff_flux = FLUX_P + FLUX_M		   
            ENDIF
            IF(XPQ(I).LT.1E-14) XPQ(I)=0.
         ENDIF
      ENDDO
      b0pom = 5.5d0 
	wghtt = DEXP(-B0POM*DABS(T2))
      RETURN
      END
