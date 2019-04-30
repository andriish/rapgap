      SUBROUTINE QCDGRIDI
	Implicit None
c test  for integration of grid values
      REAL SLO,SHI,FU,FL,TQPM,TQQB,TCCB,TBBB,TQCDC,XEPS,EPS
	Real QPMS,QQBS,CCBS,BBBS,QCDCS
      EXTERNAL FU,FL,TQPM,TQQB,TCCB,TBBB,TQCDC

      SLO = 0.0
      SHI = 1.0

      XEPS = 1e-4
      EPS = XEPS
      CALL GADAP2(SLO,SHI,FL,FU,TQPM,EPS,QPMS)
      EPS = XEPS
      CALL GADAP2(SLO,SHI,FL,FU,TQQB,EPS,QQBS)
      EPS = XEPS
      CALL GADAP2(SLO,SHI,FL,FU,TCCB,EPS,CCBS)
      EPS = XEPS
      CALL GADAP2(SLO,SHI,FL,FU,TBBB,EPS,BBBS)
      EPS = XEPS
      CALL GADAP2(SLO,SHI,FL,FU,TQCDC,EPS,QCDCS)

      write(6,*) ' QCDGRIDIN          QPM  = ',QPMS,' nb '
      write(6,*) ' QCDGRIDIN          QQB  = ',QQBS,' nb '
      write(6,*) ' QCDGRIDIN          CCB  = ',CCBS,' nb '
      write(6,*) ' QCDGRIDIN          BBB  = ',BBBS,' nb '
      write(6,*) ' QCDGRIDIN          QCDC = ',QCDCS,' nb '
      return
      end
      FUNCTION TQPM(X1,X2)
	Implicit None
      REAL TQPM,X1,X2
#include "rgqcdgrid.inc"
#include "rgpara1.inc"
#include "rgparas.inc"
#include "rgdisdif.inc"
#include "rgdiffr.inc"
#include "rglujets.inc"
#include "rgpart.inc"
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      Real Y,YD,QD,X1P,X2P,WEIG
	Integer Iy,IQ
	Double Precision ME,MP,W02,W12
	ME = PYMASS(K(1,1))
	MP = PYMASS(K(2,1))
	W02=(AM(1)+MP)**2
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))
      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
      Y = YMIN*(YMAX/YMIN)**X1
      Q2MIN=ME*ME*Y*Y/(1.D0-Y)
      IF(QMI.GT.Q2MIN) Q2MIN = QMI
      Q2MAX=Y*SSS - W12
      IF(QMA.LT.Q2MAX) Q2MAX = QMA
      Q2 = Q2MIN*(Q2MAX/Q2MIN)**X2
      IY = 0
   10 IY = IY + 1
      IF(IY+1.GE.NY) THEN
c          write(6,*) ' y  out of  grid',yy,NY,IY
         GOTO 20
      ENDIF
      IF(DBLE(Y).GT.QY(IY+1)) GOTO 10
   20 CONTINUE
      IQ = 0
   30 IQ = IQ + 1
      IF(IQ+1.GE.NQ) THEN
c          write(6,*) ' q2 out of  grid',q2,NQ,IQ
         GOTO 40
      ENDIF
      IF(Q2.GT.QQ(IQ+1)) GOTO 30
   40 CONTINUE
      YD = (DBLE(Y) - QY(IY))/(QY(IY+1) - QY(IY))
      QD = (Q2 - QQ(IQ))/(QQ(IQ+1) - QQ(IQ))
      IF(ID.EQ.1) THEN
         X1P= (QPM(IY+1,IQ)-QPM(IY,IQ))* YD +QPM(IY,IQ)
         X2P= (QPM(IY+1,IQ+1)-QPM(IY,IQ+1))* YD +QPM(IY,IQ+ 1)
      ELSEIF(ID.EQ.2) THEN
c         X1P= (QPMDF(IY+1,IQ)-QPMDF(IY,IQ))* YD +QPMDF(IY, IQ)
c         X2P= (QPMDF(IY+1,IQ+1)-QPMDF(IY,IQ+1))* YD + QPMDF(IY,IQ+1)
      ELSEIF(ID.EQ.3) THEN
c         X1P= (QPMPI(IY+1,IQ)-QPMPI(IY,IQ))* YD +QPMPI(IY, IQ)
c         X2P= (QPMPI(IY+1,IQ+1)-QPMPI(IY,IQ+1))* YD + QPMPI(IY,IQ+1)
      ENDIF
      WEIG = DBLE(Y)*DLOG(YMAX/YMIN)*Q2*DLOG(Q2MAX/Q2MIN)
      TQPM = ((X2P-X1P)*QD + X1P)*WEIG
      RETURN
      END
      FUNCTION TQQB(X1,X2)
	Implicit None
      REAL TQQB,X1,X2
#include "rgqcdgrid.inc"
#include "rgpara1.inc"
#include "rgparas.inc"
#include "rgdisdif.inc"
#include "rgdiffr.inc"
#include "rglujets.inc"
#include "rgpart.inc"
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      Real Y,YD,QD,X1P,X2P,WEIG
	Integer Iy,IQ
	Double Precision ME,MP,W02,W12
	ME = PYMASS(K(1,1))
	MP = PYMASS(K(2,1))
	W02=(AM(1)+MP)**2
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))
      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
      Y = YMIN*(YMAX/YMIN)**X1
      Q2MIN=ME*ME*Y*Y/(1.D0-Y)
      IF(QMI.GT.Q2MIN) Q2MIN = QMI
      Q2MAX=Y*SSS - W12
      IF(QMA.LT.Q2MAX) Q2MAX = QMA
      Q2 = Q2MIN*(Q2MAX/Q2MIN)**X2
      IY = 0
   10 IY = IY + 1
      IF(IY+1.GE.NY) THEN
c          write(6,*) ' y  out of  grid',yy,NY,IY
         GOTO 20
      ENDIF
      IF(DBLE(Y).GT.QY(IY+1)) GOTO 10
   20 CONTINUE
      IQ = 0
   30 IQ = IQ + 1
      IF(IQ+1.GE.NQ) THEN
c          write(6,*) ' q2 out of  grid',q2,NQ,IQ
         GOTO 40
      ENDIF
      IF(Q2.GT.QQ(IQ+1)) GOTO 30
   40 CONTINUE
      YD = (DBLE(Y) - QY(IY))/(QY(IY+1) - QY(IY))
      QD = (Q2 - QQ(IQ))/(QQ(IQ+1) - QQ(IQ))

      IF(ID.EQ.1) THEN
         X1P= (QQB(IY+1,IQ)-QQB(IY,IQ))* YD +QQB(IY,IQ)
         X2P= (QQB(IY+1,IQ+1)-QQB(IY,IQ+1))* YD +QQB(IY,IQ+ 1)
      ELSEIF(ID.EQ.2) THEN
c         X1P= (QQBDF(IY+1,IQ)-QQBDF(IY,IQ))* YD +QQBDF(IY, IQ)
c         X2P= (QQBDF(IY+1,IQ+1)-QQBDF(IY,IQ+1))* YD + QQBDF(IY,IQ+1)
      ELSEIF(ID.EQ.3) THEN
c         X1P= (QQBPI(IY+1,IQ)-QQBPI(IY,IQ))* YD +QQBPI(IY, IQ)
c         X2P= (QQBPI(IY+1,IQ+1)-QQBPI(IY,IQ+1))* YD + QQBPI(IY,IQ+1)
      ENDIF
      WEIG = DBLE(Y)*DLOG(YMAX/YMIN)*Q2*DLOG(Q2MAX/Q2MIN)
      TQQB = ((X2P-X1P)*QD + X1P)*WEIG
      RETURN
      END
      FUNCTION TCCB(X1,X2)
	Implicit None
      REAL TCCB,X1,X2
#include "rgqcdgrid.inc"
#include "rgpara1.inc"
#include "rgparas.inc"
#include "rgdisdif.inc"
#include "rgdiffr.inc"
#include "rglujets.inc"
#include "rgpart.inc"
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      Real Y,YD,QD,X1P,X2P,WEIG
	Integer Iy,IQ
	Double Precision ME,MP,W02,W12
	ME = PYMASS(K(1,1))
	MP = PYMASS(K(2,1))
	W02=(AM(1)+MP)**2
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))
      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
      Y = YMIN*(YMAX/YMIN)**X1
      Q2MIN=ME*ME*Y*Y/(1.D0-Y)
      IF(QMI.GT.Q2MIN) Q2MIN = QMI
      Q2MAX=Y*SSS - W12
      IF(QMA.LT.Q2MAX) Q2MAX = QMA
      Q2 = Q2MIN*(Q2MAX/Q2MIN)**X2
      IY = 0
   10 IY = IY + 1
      IF(IY+1.GE.NY) THEN
c          write(6,*) ' y  out of  grid',yy,NY,IY
         GOTO 20
      ENDIF
      IF(DBLE(Y).GT.QY(IY+1)) GOTO 10
   20 CONTINUE
      IQ = 0
   30 IQ = IQ + 1
      IF(IQ+1.GE.NQ) THEN
c          write(6,*) ' q2 out of  grid',q2,NQ,IQ
         GOTO 40
      ENDIF
      IF(Q2.GT.QQ(IQ+1)) GOTO 30
   40 CONTINUE
      YD = (DBLE(Y) - QY(IY))/(QY(IY+1) - QY(IY))
      QD = (Q2 - QQ(IQ))/(QQ(IQ+1) - QQ(IQ))
      IF(ID.EQ.1) THEN
         X1P= (QQBH(IY+1,IQ)-QQBH(IY,IQ))* YD +QQBH(IY,IQ)
         X2P= (QQBH(IY+1,IQ+1)-QQBH(IY,IQ+1))* YD +QQBH(IY, IQ+1)
      ELSEIF(ID.EQ.2) THEN
c         X1P= (QQBHDF(IY+1,IQ)-QQBHDF(IY,IQ))* YD + QQBHDF(IY,IQ)
c         X2P= (QQBHDF(IY+1,IQ+1)-QQBHDF(IY,IQ+1))* YD + QQBHDF(IY,IQ+1)
      ELSEIF(ID.EQ.3) THEN
c         X1P= (QQBHPI(IY+1,IQ)-QQBHPI(IY,IQ))* YD + QQBHPI(IY,IQ)
c         X2P= (QQBHPI(IY+1,IQ+1)-QQBHPI(IY,IQ+1))* YD + QQBHPI(IY,IQ+1)
      ENDIF

      WEIG = DBLE(Y)*DLOG(YMAX/YMIN)*Q2*DLOG(Q2MAX/Q2MIN)
      TCCB = ((X2P-X1P)*QD + X1P)*WEIG
      RETURN
      END
      FUNCTION TBBB(X1,X2)
	Implicit None
      REAL TBBB,X1,X2
#include "rgqcdgrid.inc"
#include "rgpara1.inc"
#include "rgparas.inc"
#include "rgdisdif.inc"
#include "rgdiffr.inc"
#include "rglujets.inc"
#include "rgpart.inc"
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      Real Y,YD,QD,X1P,X2P,WEIG
	Integer Iy,IQ
	Double Precision ME,MP,W02,W12
	ME = PYMASS(K(1,1))
	MP = PYMASS(K(2,1))
	W02=(AM(1)+MP)**2
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))
      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
      Y = YMIN*(YMAX/YMIN)**X1
      Q2MIN=ME*ME*Y*Y/(1.D0-Y)
      IF(QMI.GT.Q2MIN) Q2MIN = QMI
      Q2MAX=Y*SSS - W12
      IF(QMA.LT.Q2MAX) Q2MAX = QMA
      Q2 = Q2MIN*(Q2MAX/Q2MIN)**X2
      IY = 0
   10 IY = IY + 1
      IF(IY+1.GE.NY) THEN
c          write(6,*) ' y  out of  grid',yy,NY,IY
         GOTO 20
      ENDIF
      IF(DBLE(Y).GT.QY(IY+1)) GOTO 10
   20 CONTINUE
      IQ = 0
   30 IQ = IQ + 1
      IF(IQ+1.GE.NQ) THEN
c          write(6,*) ' q2 out of  grid',q2,NQ,IQ
         GOTO 40
      ENDIF
      IF(Q2.GT.QQ(IQ+1)) GOTO 30
   40 CONTINUE
      YD = (DBLE(Y) - QY(IY))/(QY(IY+1) - QY(IY))
      QD = (Q2 - QQ(IQ))/(QQ(IQ+1) - QQ(IQ))

      IF(ID.EQ.1) THEN
         X1P= (QQBB(IY+1,IQ)-QQBB(IY,IQ))* YD +QQBB(IY,IQ)
         X2P= (QQBB(IY+1,IQ+1)-QQBB(IY,IQ+1))* YD +QQBB(IY, IQ+1)
      ELSEIF(ID.EQ.2) THEN
c         X1P= (QQBBDF(IY+1,IQ)-QQBBDF(IY,IQ))* YD + QQBBDF(IY,IQ)
c         X2P= (QQBBDF(IY+1,IQ+1)-QQBBDF(IY,IQ+1))* YD + QQBBDF(IY,IQ+1)
      ELSEIF(ID.EQ.3) THEN
c         X1P= (QQBBPI(IY+1,IQ)-QQBBPI(IY,IQ))* YD + QQBBPI(IY,IQ)
c         X2P= (QQBBPI(IY+1,IQ+1)-QQBBPI(IY,IQ+1))* YD + QQBBPI(IY,IQ+1)
      ENDIF

      WEIG = DBLE(Y)*DLOG(YMAX/YMIN)*Q2*DLOG(Q2MAX/Q2MIN)
      TBBB = ((X2P-X1P)*QD + X1P)*WEIG
      RETURN
      END
      FUNCTION TQCDC(X1,X2)
	Implicit None
      REAL TQCDC,X1,X2
#include "rgqcdgrid.inc"
#include "rgpara1.inc"
#include "rgparas.inc"
#include "rgdisdif.inc"
#include "rgdiffr.inc"
#include "rglujets.inc"
#include "rgpart.inc"
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      Real Y,YD,QD,X1P,X2P,WEIG
	Integer Iy,IQ
	Double Precision ME,MP,W02,W12
	ME = PYMASS(K(1,1))
	MP = PYMASS(K(2,1))
	W02=(AM(1)+MP)**2
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))
      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
      Y = YMIN*(YMAX/YMIN)**X1
      Q2MIN=ME*ME*Y*Y/(1.D0-Y)
      IF(QMI.GT.Q2MIN) Q2MIN = QMI
      Q2MAX=Y*SSS - W12
      IF(QMA.LT.Q2MAX) Q2MAX = QMA
	
      Q2 = Q2MIN*(Q2MAX/Q2MIN)**X2
      IY = 0
   10 IY = IY + 1
      IF(IY+1.GE.NY) THEN
c          write(6,*) ' y  out of  grid',yy,NY,IY
         GOTO 20
      ENDIF
      IF(DBLE(Y).GT.QY(IY+1)) GOTO 10
   20 CONTINUE
      IQ = 0
   30 IQ = IQ + 1
      IF(IQ+1.GE.NQ) THEN
c          write(6,*) ' q2 out of  grid',q2,NQ,IQ
         GOTO 40
      ENDIF
      IF(Q2.GT.QQ(IQ+1)) GOTO 30
   40 CONTINUE
      YD = (DBLE(Y) - QY(IY))/(QY(IY+1) - QY(IY))
      QD = (Q2 - QQ(IQ))/(QQ(IQ+1) - QQ(IQ))


      IF(ID.EQ.1) THEN
         X1P= (QCDC(IY+1,IQ)-QCDC(IY,IQ))* YD +QCDC(IY,IQ)
         X2P= (QCDC(IY+1,IQ+1)-QCDC(IY,IQ+1))* YD +QCDC(IY, IQ+1)
      ELSEIF(ID.EQ.2) THEN
c         X1P= (QCDCDF(IY+1,IQ)-QCDCDF(IY,IQ))* YD + QCDCDF(IY,IQ)
c         X2P= (QCDCDF(IY+1,IQ+1)-QCDCDF(IY,IQ+1))* YD + QCDCDF(IY,IQ+1)
      ELSEIF(ID.EQ.3) THEN
c         X1P= (QCDCPI(IY+1,IQ)-QCDCPI(IY,IQ))* YD + QCDCPI(IY,IQ)
c         X2P= (QCDCPI(IY+1,IQ+1)-QCDCPI(IY,IQ+1))* YD + QCDCPI(IY,IQ+1)
      ENDIF

      WEIG = DBLE(Y)*DLOG(YMAX/YMIN)*Q2*DLOG(Q2MAX/Q2MIN)

      TQCDC = ((X2P-X1P)*QD + X1P)*WEIG
      RETURN
      END
