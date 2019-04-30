*###################################################################
* H1 QCD/Regge Fits
*
* Note on xpq: elements 0,1,2,3,4,5,6 are g,d,u,s,c,b,t
*              ie xpq(1)=d !
* Options
* -------
*
* NPOM = -10 H1 FINAL 1994 FIT (Flux Part)
* NG   = -10 H1 FINAL 1994 FIT (Fit 1, Quarks Only)
* NG   = -11 H1 FINAL 1994 FIT (Fit 2, Quarks + Flat Gluon)
* NG   = -12 H1 FINAL 1994 FIT (Fit 3, Quarks + Flat Gluon)
*####################################################################

      SUBROUTINE H1QCDFIT(BETA_IN,Q2_IN,XPQ,X_POM_IN,T2_IN)

      IMPLICIT REAL*8 (A-G,O-Z)
      REAL*4 T2_IN,X_POM_IN,BETA_IN,Q2_IN,XPQ(-6:6)
      DIMENSION XPQP(-6:6),XPQM(-6:6)
      LOGICAL FIRST
      INTEGER ICALL
      DATA ICALL/0/
      LOGICAL FLAG
      DATA FLAG/.FALSE./
#include "rgdiffr.inc"
#include "rgnlodiff.inc"
#include "rgpdiss.inc"

C Function for integral over t
      T_INT(A_TMIN,A_TMAX,B)=(-DEXP(-B*A_TMAX)/B) + (DEXP(-B*A_TMIN)/B)
      DATA FIRST/.TRUE./
      Icall = icall + 1
      if(icall.eq.1) then
         first=.true.
      else
         first=.false.
      endif
C      IF(FIRST) WRITE(6,*)'NPOM',NPOM
C Set Parameters for fits
      IF(NPOM.GT.-10) THEN
*! alpha_prime for the pomeron
         APPOM=0.0D0
*! b0 for pomeron
         B0POM=6.0D0
*! alpha_prime for the meson
         APMES=0.86D0
*! b0 for meson
         B0MES=1.4D0
      ELSEIF(NPOM.LE.-10 .AND. NPOM.GE.-12) THEN
         APPOM=0.26D0
         B0POM=4.6D0
         APMES=0.9D0
         B0MES=2.0D0
      ELSE
         WRITE(6,*)'H1QCDFIT: Unknown fit',NPOM
         STOP
      ENDIF
c for pdiss set fixed B slope
      IF(NFRAG.EQ.20) then
         APPOM=0.0D0
         B0POM=4.0D0
         APMES=0.0D0
         B0MES=0.0D0
	endif
      IF(NPOM.LE.-10 .AND. NPOM.GE.-12) THEN
         CP=1.D0
*! Note: 2*alpha(0)-1 -> this is alpha(0)=1.20
         PP=1.40D0
*!  "        "             "   " alpha(0)=0.57
         PM=0.14D0
         IF(NG.EQ.-10) THEN
            CM=0.011816D0
         ELSEIF(NG.EQ.-11) THEN
            CM=0.0086814D0
         ELSEIF(NG.EQ.-12) THEN
            CM=0.00863D0
         ELSEIF(NG.EQ.-13) THEN
            CM=0.011816D0
         ELSEIF(NG.EQ.-14) THEN
            CM=0.0086814D0
         ELSEIF(NG.EQ.-15) THEN
            CM=0.00863D0
         ELSE
            WRITE(6,*)'H1QCDFIT - Unknown 1994 Final fit',ABS(NG)-9
            STOP
         ENDIF
      ENDIF

C Inform user of folly
      IF(FIRST) THEN
         WRITE(6,*)'#############################################'
         WRITE(6,*)'#           H1QCD fit Selected              #'
         WRITE(6,*)'#       for Q2 > 75  fix pdf at Q2 = 75     #'
         WRITE(6,*)'#       for Q2 < 4.5 fix pdf at Q2 = 4.5    #'
         WRITE(6,*)'#   for beta < 0.04  fix pdf at beta = 0.04 #'
         IF(NPOM.LE.-10 .AND. NPOM.GE.-12) THEN
            WRITE(6,*)'# Final 1994 Fits                           #'
            WRITE(6,*)'# FIT number in paper is ',ABS(NG)-9
            IF(NPOM.EQ.-10) WRITE(6,*)'# Pomeron Part Only'
            IF(NPOM.EQ.-11) WRITE(6,*)'# Meson Part Only'
            IF(NPOM.EQ.-12) WRITE(6,*)'# Pomeron + Meson'
         ELSE
            WRITE(6,*)'# NPOM=',NPOM,'not implemented              #'
            STOP
         ENDIF
         WRITE(6,*)'#############################################'
      ENDIF

C Input quantities are REAL*4 -> Convert to REAL*8
      BETA =DBLE(BETA_IN)
      if(BETA.LT.0.04) BETA = 0.04D0
      Q2   =DBLE(Q2_IN)
      if(Q2.GT.75d0) Q2=75.d0
      if(Q2.LT.4.5d0) Q2=4.5d0
      X_POM=DBLE(X_POM_IN)
      T2   =DBLE(T2_IN)
      IF(NPOM.NE.-3) THEN
C Initialise PDF Pion
         IF(NPOM.EQ.-10 .OR. NPOM.EQ.-11 .OR. NPOM.EQ.-12) THEN
*! LO-Owens Set 1
c        write(6,*) ' LHApdf called '
          CALL PION_INIT(0.D0,0.D0,211.D0) ! LO-Owens Set 1
c        write(6,*) ' pdflib called '
c          CALL PION_INIT(2.D0,1.D0,1.D0) ! LO-Owens Set 1
         ENDIF
      ENDIF

C Final 1994 Fits from paper
      IF(NPOM.LE.-10 .AND. NPOM.GE.-12) THEN
         IF(FIRST) THEN
            IFIT=ABS(NG)-9
         ELSE
            IFIT=0
         ENDIF
         CALL QCD_1994(BETA,Q2,XPQP,IFIT)
c use Q instead of Q2 for PDFLIB 6.4.98
ccc         CALL PFTOPDG( BETA,DSQRT(Q2),XPQM)
         CALL PFTOPDG( BETA,Q2,XPQM)
	   do i=1,6
	     xpqm(i) =  (xpqm(-i)+xpqm(i))/2.
	     xpqm(-i) = xpqm(i)
	   enddo
         XPR=0.003D0
         ATMIN = (0.93827231D0*XPR)**2/(1.D0-XPR)
         DP = T_INT(ATMIN,1.D0,B0POM+2.D0*APPOM*LOG(1.D0/XPR)) *XPR**(-
     +   PP)
         DM = T_INT(ATMIN,1.D0,B0MES+2.D0*APMES*LOG(1.D0/XPR)) *XPR**(-
     +   PM)
         FLUX_P = 0.D0
         FLUX_M = 0.D0
         CALL FLUXH1(X_POM,T2,APPOM,B0POM,CP,PP,FLUX_P)
         CALL FLUXH1(X_POM,T2,APMES,B0MES,CM,PM,FLUX_M)
         FLUX_P = FLUX_P/(DP*XPR)
         FLUX_M = FLUX_M/(DM*XPR)
      ENDIF

      DO I=-6,6
         XPQ(I)=0.
         IF(ABS(I).LT.6) THEN
            IF(NPOM.EQ.-10 .OR. NPOM.EQ.-3) THEN
               XPQ(I)=REAL(XPQP(I)*FLUX_P)
		   Diff_flux = FLUX_P
            ELSEIF(NPOM.EQ.-11 ) THEN
               XPQ(I)=REAL(XPQM(I)*FLUX_M)
		   Diff_flux = FLUX_M		   
            ELSEIF(NPOM.EQ.-12) THEN
               XPQ(I)=REAL(XPQM(I)*FLUX_M+XPQP(I)*FLUX_P)
		   Diff_flux = FLUX_P + FLUX_M		   
            ENDIF
		
            IF(XPQ(I).LT.1E-10) XPQ(I)=0.
         ENDIF
         IF(FLAG) WRITE(6,*)'I,XPQ(I)',I,XPQ(I)
      ENDDO
	wghtt = DEXP(-B0POM*DABS(T2))
c      write(6,*) wghtt,b0pom,t2
      FIRST=.FALSE.

      RETURN
      END

      subroutine pftopdg(dx,dscale,dxpdf)
      double precision dx,dscale,dupv,ddnv,dusea,ddsea,
     &  dstr,dchm,dbot,dtop,dgl,dxpdf(-6:6)
      call structm(dx,dscale,dupv,ddnv,dusea,ddsea,dstr,dchm,dbot,dtop,dgl)
c Convert flavour convention of PDFLIB to PDG convention
      dxpdf(0) = dgl
      dxpdf(1) = ddnv + ddsea
      dxpdf(2) = dupv + dusea
      dxpdf(3) = dstr
      dxpdf(4) = dchm
      dxpdf(5) = dbot
      dxpdf(6) = dtop
      dxpdf(-1) = ddsea
      dxpdf(-2) = dusea
      dxpdf(-3) = dstr
      dxpdf(-4) = dchm
      dxpdf(-5) = dbot
      dxpdf(-6) = dtop
      return
      end
      
      
