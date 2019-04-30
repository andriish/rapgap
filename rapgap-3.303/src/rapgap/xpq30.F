      Subroutine XPQ30(X,Q2,XPQ)
      Implicit None
      Integer I
      REAL X,Q2,XPQ
      Double Precision XD,Q2D,SIGOUT,ESIGOUT,WEIOUT
      Double Precision pdfnorm
      Dimension XPQ(-6:6)
#include "rgparam.inc"
      DOUBLE PRECISION POM,PIM
      COMMON/WEIGHT1/ POM,PIM
      Double Precision GEV2NB
      DATA GEV2NB/.3893857D+3/

      XD = X
      Q2D = Q2
      Do I=-6,6
         XPQ(I) = 0.
      Enddo

      Call SATSIGTOT(Q2D,XD,SIGOUT,ESIGOUT,WEIOUT)
c     write(6,*) 'XPQ30:Q2D,XD,SIGOUT,ESIGOUT,WEIOUT',
c     &  Q2D,XD,SIGOUT,ESIGOUT,WEIOUT
c factor 9 comes from 1/eq**2 for d_quark with eq=1/3
ccc   XPQ(1) = SIGOUT*Q2D*9.D0/4.D0/PI/PI/ALPH
      pdfnorm=1./2.d0/9.d0
      pdfnorm=Q2D*9.D0/12.D0/4.D0/PI/PI/ALPH/GEV2NB
      XPQ(1) = SIGOUT*pdfnorm
      XPQ(-1) = SIGOUT*pdfnorm
      XPQ(2) = SIGOUT*pdfnorm
      XPQ(-2) = SIGOUT*pdfnorm
      XPQ(3) = SIGOUT*pdfnorm
      XPQ(-3) = SIGOUT*pdfnorm
      POM = WEIOUT
c      write(6,*) 'XPQ30: XPQ(1) ',XPQ(1),' weight ',weiout,' x,q2 ',x,q2

      RETURN
      END
