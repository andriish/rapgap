      Subroutine RGSATRAP(X,WDIF)
      Implicit NONE
#include "rgpart.inc"
#include "rgpara1.inc"
#include "rgdiffr.inc"
      Double Precision X(20)
      Double Precision WDIF

      Double Precision  Stot, Q2s, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2s, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision    DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT
      Double Precision    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT
      COMMON /CCOMPCROS/  DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT,
     &                    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT

      Integer         NQUAD
      COMMON /CNQUAD/ NQUAD

      Double Precision  EPS
      Common  /mysatrap/EPS
      Double Precision W02,W12
      DOUBLE PRECISION ME,MP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Get Parameters of the Wuesthoff Model

      Call SATURPAR
C -----------------------------------------------------

      WDIF = 0.
      ME =0.511e-3
      MP =0.938


      Stot = SSS

c calculate kinematically alowed y range:
      W02=(1.+MP)**2
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))


      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
c     write(6,*) 'rgsatrap ',ymin,ymax

c Mx limits moved to    RGQ2W2Mx

      NQUAD = 24

c     write(6,*) ' in RGSATRAP '

c Compute diffractive sigma's by MC computation

      Call RGQ2W2Mx(X)
      IF(AJAC.LE.0) RETURN
      Call INPUTKINE

      Call COMPCROSS

      EPS = 2.*(1.-YBj)/(1.+(1.-YBj)**2)

      DLTOT = DLQQT+DLQQGT + EPS*(DLQQL)
      DLCTOT= DCQQT+DCQQGT + EPS*(DCQQL+DCQQGL) + DLTOT

      IF(DLCTOT.LE.0.0) DLCTOT=0.
      WDIF = DLCTOT

      Return
      End
