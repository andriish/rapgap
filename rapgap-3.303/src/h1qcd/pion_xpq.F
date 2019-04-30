C Construct PDFs for pion or other meson
      SUBROUTINE PION_XPQ(XX,Q2,XPQ)
      IMPLICIT REAL*8 (A-G,O-Z)
      DIMENSION XPQ(-6:6)
      DOUBLE PRECISION XX,QQ,UPV,DNV,SEA,STR,CHM,BOT,TOP,GLU

C PDFlib only contains data for the negative pion, so here the dval and
C uval are swapped to give a positive pion
*! ud isospin 0 meson
      ISF=2
C      ISF=1 ! pi+

      QQ = DSQRT(Q2)

C Initialise
      DO I=-6,6
         XPQ(I) = 0.D0
      ENDDO

      IF(XX.GE.1. .OR. XX.LE.0.) RETURN
c       write(6,*) 'pion_xpq: X,Q', xx,QQ
c      CALL STRUCTF(XX,QQ,UPV,DNV,SEA,STR,CHM,BOT,TOP,GLU)
      CALL STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)

C Positive Pion
      IF(ISF.EQ.1) THEN
*! Gluon
         XPQ(0)= GLU
*! D
         XPQ(1)= UPV+USEA
*! Dbar
         XPQ(-1)= USEA
*! U
         XPQ(2)= DNV+DSEA
*! Ubar
         XPQ(-2)= DSEA
*! S
         XPQ(3)= STR
*! Sbar
         XPQ(-3)= STR
*! C
         XPQ(4)= CHM
*! Cbar
         XPQ(-4)= CHM
*! B
         XPQ(5)= BOT
*! Bbar
         XPQ(-5)= BOT
*! T
         XPQ(6)= TOP
*! Tbar
         XPQ(-6)= TOP
      ELSEIF(ISF.EQ.2) THEN
C Flavour singlet meson (f2)
         VAL = UPV/2.D0
*! Gluon
         XPQ(0)= GLU
*! D
         XPQ(1)= VAL+USEA
*! Dbar
         XPQ(-1)= VAL+USEA
*! U
         XPQ(2)= VAL+USEA
*! Ubar
         XPQ(-2)= VAL+USEA
*! S
         XPQ(3)= STR
*! Sbar
         XPQ(-3)= STR
*! C
         XPQ(4)= CHM
*! Cbar
         XPQ(-4)= CHM
*! B
         XPQ(5)= BOT
*! Bbar
         XPQ(-5)= BOT
*! T
         XPQ(6)= TOP
*! Tbar
         XPQ(-6)= TOP
      ENDIF

      RETURN
      END
