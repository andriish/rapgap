      FUNCTION D_XGX(XX,QQ2)
	Implicit None
      DOUBLE PRECISION XGX,X,XX,QQ2,DELTA,DXGX,RERR,D_XGX
      DOUBLE PRECISION DQ2
      COMMON/gluon1/X
      EXTERNAL XGX
      DELTA = 1.D-1*QQ2
      RERR = -9999.D0
      DXGX = -99999.D0
      DQ2 = QQ2
      X = XX
      CALL DFRIDR(XGX,DQ2,DELTA,DXGX,RERR)
      IF(RERR.EQ.-9999.D0) THEN
         write(6,*) ' D_XGX: error in calculating derivative of xg(x,'
     +   //'q2)'
         write(6,*) ' error = ',RERR
         DXGX = -99999.D0
      ENDIF
      IF(RERR/dxgx.GT.5.D-2) THEN
         write(6,*) ' D_XGX: large error for derivative of xg(x,'
     +   //'q2)'
         write(6,*) ' error = ',RERR/dxgx,' larger than 0.01 '
         write(6,*) ' with RERR = ',RERR,' and dxgx = ',dxgx
	   write(6,*) ' results are not reliable'
	   write(6,*) ' check selected set of pdf"s'
	   write(6,*) ' are they valid for Q2 = ',DQ2,' and x = ',XX
	   write(6,*) ' stop the program; derivative set to 0.0000 '
         DXGX = 0.D0
      ENDIF
c      write(6,*) ' D_XGX: DFRIDR DXGX,DELTA,RERR ',DXGX,DELTA,RERR
      D_XGX = DXGX
      RETURN
      END


      FUNCTION XGX(QQ2)
	Implicit none
      DOUBLE PRECISION X,XGX,QQ2
      COMMON/gluon1/X
      DOUBLE PRECISION XXX,Q2,XGLU
	EXTERNAL XGLU
      Q2=QQ2
      XXX = X
      XGX = XGLU(XXX,QQ2)
      RETURN
      END
      FUNCTION XGLU(XX,QQ2)
	Implicit none
      DOUBLE PRECISION XX,QQ2,XGLU
      DOUBLE PRECISION XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     +VAL(20)
      LOGICAL PDFFIRST
      COMMON /W50516/PDFFIRST
      CHARACTER*20 PARM(20)
#include "rgpypars.inc"
      Integer MST51,NPDF
      DATA NPDF/0/,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU/9*0.D0/
      QQ=SQRT(MAX(0.D0,QQ2))
      XXX = XX
      PARM(1) = 'NPTYPE'
      MST51 = MSTP(51)
      VAL(1) = INT(MST51/1000000)

      PARM(2) = 'NGROUP'
      VAL(2) = DBLE(MOD(MST51,1000000)/1000)
      PARM(3) = 'NSET'
      VAL(3) = DBLE(MOD(MOD(MST51,1000000),1000))
CNEW
      NPDF=NPDF+1
      PDFFIRST=.FALSE.
      IF(NPDF.LE.1) then
        PDFFIRST=.TRUE.
      Endif
ccc for LHAPDF 5.7            
      PARM(1) = 'DEFAULT'
      VAL(1) = mstp(51)
ccc

c call PDFSET each time, because when DIF with pi structure function is
c also selected one would get in confusion....
      CALL PDFSET(PARM,VAL)
      CALL STRUCTM(XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
      XGLU = GLU
      RETURN
      END
