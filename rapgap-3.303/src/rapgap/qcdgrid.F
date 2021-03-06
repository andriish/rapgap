      SUBROUTINE QCDGRID
      IMPLICIT None
      Double Precision X,xv,wmax
      DIMENSION X(20)
      Integer NDIM,NPOIN,NDIMEN,IMIX
      COMMON/DIVO/ NDIM,NPOIN
      COMMON /XVAL/ XV(20),NDIMEN
      REAL RQPM,RQQB,RQQBC,RQQBB,RQCDC,FL,FU,PQCD,XEPS,EPS,PQCDI
      COMMON /OALPHAS/ WMAX,IMIX
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
#include "rgqcdgrid.inc"
#include "rgraher.inc"
#include "rgpart.inc"
#include "rghs45.inc"
#include "rgrapgki.inc"
#include "rgpara.inc"
#include "rgpara1.inc"
#include "rgdisdif.inc"
#include "rgluco.inc"
#include "rgdiffr.inc"
      Double Precision PYMASS
      REAL SLO,SHI
      EXTERNAL FU,FL,PQCD,PQCDI,PYMASS
      double precision y,xxx,RSAFETY,xpo_max,xpo_min
      double precision RSUM,PQPM
      Integer idiro,iproo,iheraco,nflavo,ndimeno,ndimo,ngo,npomo
      Integer nny,nnq,nnpom,ipom,IHERPYSO,iy,iq,i,iprint
      DATA IPRINT/0/
      XEPS = 0.001
      XEPS = 1e-4
      SLO = 0.0
      SHI = 1.0
c      write(6,*) ' rapgap check weights ',IQCDGRID
c      write(6,*) idir,iherpys
      IDIRO = IDIR
      IPROO = IPRO
      IHERACO = IHERAC
      NFLAVO =NFLAV
      NDIMENO =NDIMEN
      NDIMO = NDIM
      IHERPYSO = IHERPYS
      NGO = NG
      NPOMO = NPOM
      IDIR = 1
	IDIRINT = IDIR
      IHERAC = 0
      NDIMEN = 2
      NDIM =2
      IF(IDIRO.EQ.1) THEN
         ID = 1
      ELSEIF(IDIRO.EQ.0.AND.(NPOM.EQ.20.OR.NPOM.EQ.21)) THEN
         ID = 3
      ELSE
         ID = 2
      ENDIF
C select leading order or higher order process
      IF(IDISDIF.GE.3.AND.IDIRO.EQ.1) THEN
         write(6,*) ' this option is  not impelmented --> STOP '
         STOP
      ENDIF
      IF(IDISDIF.GE.1.AND.IDIRO.EQ.1) THEN
         IHERPYS = 0
         ID = 1
      ENDIF
      write(6,*) ' Now calcualting order alpha_s matrix elements '
      write(6,*) ' Please be patient, for diffraction it can take long'
      IF(IPRINT.EQ.1) then
         IF(ID.EQ.1) THEN
            WRITE(6,10000)
10000 FORMAT('|   y     |   q2    |  qpm    |  qqb    |  qcdc   |',
     + '  ccb    |   bbb   |')
         ELSE
            WRITE(6,10100)
10100 FORMAT('|   y     |   q2    |  x_pom  |  qpm    |  qqb    |',
     + '  qcdc   |  ccb    |   bbb   |')
         ENDIF
      endif
c choose y and q2
   10 CONTINUE
      NNY = NY - 1
      NNQ = NQ - 1
      NNPOM = NPPO -1
      IF(ID.EQ.1) THEN
         NNPOM=0
      ELSE
         XPO_MAX=1.D0-XF
         XPO_MIN=Q2MIN/YMAX/SSS
c        write(6,*) T2GKI
         T2GKI = -0.1
         XFGKI = 0.01
      ENDIF
c     write(6,*) ' qcdgrid NNpom,xpo_min,xpo_max ',NNpom,xpo_min,xpo_max
      DO 70 Ipom = 0,NNPOM
         DO 60 IY = 0,NNY
            DO 50 IQ = 0,NNQ
               DO 20 I=1,10
   20          X(I) = 0.D0
               X(1) = DFLOAT(IY)/DFLOAT(NNY)
               X(2) = DFLOAT(IQ)/DFLOAT(NNQ)
               X(3) = DFLOAT(IPOM)/DFLOAT(NNPOM)
               DO 30 I=1,20
   30          XV(I) = X(I)
               IDIR = 1
		   IDIRINT = IDIR
               IHERAC = 0
               NDIMEN = 2
               IMIX = 0
               IWEI = 0
               IPRO = 12
               KPA = 1
ccccc            NFLAV = 3

               Y = YMIN*(YMAX/YMIN)**X(1)
               YY = SNGL(y)
               Q2 = Q2MIN*(Q2MAX/Q2MIN)**X(2)
               IF(ID.EQ.1) THEN
                  XPOM = 1.D0
               ELSE
                  XPOM = XPO_MIN*(XPO_MAX/XPO_MIN)**X(3)
               ENDIF
cc       write(6,*) x(1),y,ymin,ymax,x(2),q2,q2min,q2max
               xxx = q2/y/sss
               IF(xxx.ge.1.) THEN
c                 write(6,*) 'iY,Iq',iy,iq
                  Q2=y*sss*0.99D0
                  goto 50
               endif
               RQPM = SNGL(PQPM(Y,Q2))
               IMIX = 1
               IWEI = 1
               KPA = 1
               RQQB=0.0
               RQQBC=0.0
               RQCDC=0.0
               if(RQPM.LE.0.0) THEN
                  RQPM = 1.e-10
                  GOTO 40
               ENDIF
               IF(IHF.EQ.0) THEN
                  IPRO = 13
c calculate maximum of diff x section
C........ gamma gluon fusion
                  WMAX = 0.D0
                  AM(1) = PYMASS(1)
                  AM(2) = AM(1)
                  EPS = XEPS
                  IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
                     CALL GADAP2(SLO,SHI,FL,FU,PQCD,EPS,RQQB)
                  else
                     CALL GADAP(SLO,SHI,PQCDI,EPS,RQQB)
                  ENDIF
c           write(6,*) ' gamma ghluon fusion'
c           CALL DULIST(1)
C........ QCDC
                  IPRO = 15
c calculate maximum of diff x section
                  WMAX = 0.D0
                  AM(1) = PYMASS(1)
                  AM(2) = AM(1)
                  EPS = XEPS
c           write(6,*) 'event: IPRO =',IPRO,KPA
                  IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
                     CALL GADAP2(SLO,SHI,FL,FU,PQCD,EPS,RQCDC)
                  else
                     CALL GADAP(SLO,SHI,PQCDI,EPS,RQCDC)
                  endif
c           write(6,*) ' QCD compton'
c           CALL DULIST(1)
               ENDIF
C........ gamma gluon fusion charm
               IPRO = 14
               KPA = 4
cccc            IF(NFLAVO.GE.4.AND.IHF.EQ.0) THEN
               IF(NFLAVO.GE.4) THEN
C........ gamma gluon fusion charm
                  AM(1) = PYMASS(KPA)
                  AM(2) = AM(1)
c calculate maximum of diff x section
                  WMAX = 0.D0
                  EPS = XEPS
c            write(6,*) 'qcdgrid: IPRO =',IPRO,KPA,am(1)
                  IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
                     CALL GADAP2(SLO,SHI,FL,FU,PQCD,EPS,RQQBC)
                  else
                     CALL GADAP(SLO,SHI,PQCDI,EPS,RQQBC)
                  endif
                  IF(IHF.EQ.1.AND.IHFLA.EQ.5) RQQBC=0
                  RQQBC=RQQBC
                  IF(NFLAVO.GE.5) THEN
c               write(6,*) 'qcdgrid bbar',nflavo
                     KPA = 5
                     AM(1) = PYMASS(5)
                     AM(2) = AM(1)
                     EPS = XEPS
                     IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
                        CALL GADAP2(SLO,SHI,FL,FU,PQCD,EPS,RQQBB)
                     else
                        CALL GADAP(SLO,SHI,PQCDI,EPS,RQQBB)
                     endif
                     IF(IHF.EQ.1.AND.IHFLA.EQ.4) RQQBB=0
                  ENDIF
               ENDIF
   40          CONTINUE
c            write(6,*) 'ccbar RQPM  RQQBC',RQPM,RQQBC
c           write(6,*) ' gamma gluon fusion c cbar'
c           CALL DULIST(1)
               QY(IY+1) = DBLE(YY)
               QQ(IQ+1) = Q2
               QPOM(IPOM+1) = XPOM
c            write(6,*) ' y = ',yy,' Q2 = ',Q2
c            write(6,*) ' RQPM = ',RQPM,' RQQB = ',RQQB,'RQCDC = ',RQCDC,
c     +      'RQQBC = ',RQQBC
c            write(6,*) 'QCDGRID',IHERPYS
               IF(ID.EQ.1) THEN
                  QPM(IY+1,IQ+1) = DBLE(RQPM)
                  QQB(IY+1,IQ+1) = DBLE(RQQB)
                  QCDC(IY+1,IQ+1) = DBLE(RQCDC)
                  QQBH(IY+1,IQ+1) = DBLE(RQQBC)
                  QQBB(IY+1,IQ+1) = DBLE(RQQBB)
               ELSEIF(ID.EQ.2) THEN
                  QPMDF(IY+1,IQ+1,IPOM+1) = DBLE(RQPM)
                  QQBDF(IY+1,IQ+1,IPOM+1) = DBLE(RQQB)
                  QCDCDF(IY+1,IQ+1,IPOM+1) = DBLE(RQCDC)
                  QQBHDF(IY+1,IQ+1,IPOM+1) = DBLE(RQQBC)
                  QQBBDF(IY+1,IQ+1,IPOM+1) = DBLE(RQQBB)
               ELSEIF(ID.EQ.3) THEN
                  QPMPI(IY+1,IQ+1,IPOM+1) = DBLE(RQPM)
                  QQBPI(IY+1,IQ+1,IPOM+1) = DBLE(RQQB)
                  QCDCPI(IY+1,IQ+1,IPOM+1) = DBLE(RQCDC)
                  QQBHPI(IY+1,IQ+1,IPOM+1) = DBLE(RQQBC)
                  QQBBPI(IY+1,IQ+1,IPOM+1) = DBLE(RQQBB)
               ENDIF
               RSAFETY = DBLE(0.05*RQPM)
               RSUM = DBLE(RQQB+RQCDC+RQQBC+RQQBB)
               if(IPRINT.EQ.1.and.RQPM.gt.2.e-10) THEN
c               if(IPRINT.EQ.1) THEN
                  IF(ID.EQ.1) THEN
                     WRITE(6,10200) yy,Q2,RQPM,RQQB,RQCDC,RQQBC,RQQBB
                  ELSE
                     WRITE(6,10300) yy,Q2,XPOM,RQPM,RQQB,RQCDC,RQQBC,
     +               RQQBB
                  ENDIF
               ENDIF
               IF(RSUM.GT.(DBLE(RQPM)+RSAFETY).and.IHF.EQ.0) THEN
                  IF(ID.EQ.1) THEN
                     WRITE(6,10000)
                     WRITE(6,10200) yy,Q2,RQPM,RQQB,RQCDC,RQQBC,RQQBB
                  ELSE
                     WRITE(6,10100)
                     WRITE(6,10300) yy,Q2,XPOM,RQPM,RQQB,RQCDC,RQQBC,
     +               RQQBB
                  ENDIF
                  write(6,*) ' QCDGRID: RSUM > RQPM ',RSUM,RQPM
                  write(6,*) ' including 5% safety margin '
                  write(6,*) ' increase pt cutoff PT2CUT '
                  write(6,*) ' program stopped now '
                  STOP
               ENDIF
10200 FORMAT(1X,7(E8.3,' |'))
10300 FORMAT(1X,8(E8.3,' |'))
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
c      call qcdgridi
      IF(IDISDIF.GE.1.AND.ID.EQ.1) THEN
         ID = 2
         IHERPYS = 1
         write(6,*) ' now  calculate inclusive DIF ',IHERPYSO
         GOTO 10
      ELSEIF(IDISDIF.EQ.2.AND.ID.EQ.2) THEN
         ID = 3
         IHERPYS = 1
         write(6,*) ' now  calculate inclusive PI ',IHERPYSO
         NG = 20
         NPOM = 20
         GOTO 10
      ENDIF
      IDIR = IDIRO
	IDIRINT = IDIR
      IPRO = IPROO
      IHERAC = IHERACO
      NFLAV = NFLAVO
      NDIMEN = NDIMENO
      NDIM = NDIMO
      NG = NGO
      NPOM = NPOMO
      IMIX = 0
      IWEI = 0
      IHERPYS = IHERPYSO
      RETURN
      END
