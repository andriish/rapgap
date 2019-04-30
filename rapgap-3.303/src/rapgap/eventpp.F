      SUBROUTINE EVENTPP
      IMPLICIT NONE
      Integer NDIM,NPOIN
      COMMON/DIVO/ NDIM,NPOIN
#include "rgeffic.inc"
#include "rgpara.inc"
#include "rgpara1.inc"
#include "rgluco.inc"
#include "rgpart.inc"
#include "rgludat1.inc"
#include "rgrapgki.inc"
#include "rgdiffr.inc"
#include "rgfull.inc"
#include "rgdisdif.inc"
#include "rghs45.inc"
#include "rgraher.inc"
#include "rgrgridf2.inc"
#include "rglujets.inc"
#include "rgherup.inc"
#include "rgcolcon.inc"
      Double Precision  Stot, Q2_henry, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2_henry, W2, Mx2, AJAC
      Double Precision    KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT
      COMMON /C3BODVAR/   KtG2, QtG2, SM2, ZW, KW2, WT3, YCUT

      Double Precision X(20)
      Double Precision XG(20)
      Integer LST,IRES,MSTA
      COMMON /EPPARA/ LST(30),IRES(2)
      REAL PARA
      COMMON /ARDAT1/ PARA(40),MSTA(40)
      Double Precision WMAX
      Integer IMIX
      COMMON /OALPHAS/ WMAX,IMIX

      Double Precision WPA
      Integer NQPM,NQQB,NQQBC,NQQBB,NQCDC
      COMMON /QCDEV/ NQPM,NQQB,NQQBC,NQQBB,NQCDC
      Double Precision XV
      Integer NDIMEN
      COMMON /XVAL/ XV(20),NDIMEN
      Double Precision draprn
      DOUBLE PRECISION MP
      DOUBLE PRECISION POMAX,PIMAX,WTMA
      COMMON/WEIGHT/ POMAX(NBX,NBQ2),PIMAX(NBX,NBQ2),WTMA
      Integer IGENFL,IVM
      COMMON/GENWEI/IGENFL
      COMMON/VMESON/IVM
      Integer NDIS,NQPMDS,NQQBDS,NQQBCDS,NQQBBDS,NQCDCDS,
     +              NDIF,NQPMDF,NQQBDF,NQQBCDF,NQQBBDF,NQCDCDF,
     +              NPI,NQPMPI,NQQBPI,NQQBCPI,NQQBBPI,NQCDCPI,
     +              NRPA,NRPB,NRPC,NRPD,NRPE,NRPF,NRPG,NRPH,NPRT
      COMMON/NEVOUT/NDIS,NQPMDS,NQQBDS,NQQBCDS,NQQBBDS,NQCDCDS,
     +              NDIF,NQPMDF,NQQBDF,NQQBCDF,NQQBBDF,NQCDCDF,
     +              NPI,NQPMPI,NQQBPI,NQQBCPI,NQQBBPI,NQCDCPI,
     +              NRPA,NRPB,NRPC,NRPD,NRPE,NRPF,NRPG,NRPH,NPRT
      Double Precision PYPAR,PYVAR
      Integer IPY
      COMMON /RGPYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      Integer LERR
      COMMON/ERR/LERR(100)
      REAL XHF,Q2HF
      Integer KPHF
      COMMON /HEAVYF/XHF,Q2HF,KPHF
      Integer IDEBUG
      COMMON/ INTERN/IDEBUG
      Double Precision    remPARJ32
      LOGICAL FIRST

      Integer IQPML1,IQQBL1,IQQBCL1,IQCDCL1,IQPML0
      Integer IQQBL0,IQCDCL0,IQPMH1,IQQBH1,IQQBCH1
      Integer IQCDCH1,IQPMH0,IQQBH0,IQCDCH0
      Integer NPRIN,IDEBUGO,IPROO,IDIRO,NGO,NPOMO,KPAO,KINT22
      Integer IHERPYSO,IHFO
      Integer I
      Double Precision WDUM,SMALL,T2MAXO,XBJ,XRMIN,XR
      REAL XPQ(-6:6),XPQ2(-6:6),XPQ2_DIFF(-6:6),XPRT,Q2T

      Double Precision RNTEST
      Double Precision T2MIN,XP2T,T2,T2MN,WTDIST,WTT,WTRN
      Double Precision FXN1,DCORR,XD,QD,X1P,X2P,RF2DIS,RF2DIF,RF2PI
      Double Precision Fact
      
      Integer nflavt
      Integer NTRY,IX,IQ,IERR,IKEEP
      Integer NTRMAX
      DATA NTRMAX/1000/
      
      EXTERNAL draprn
      Integer nreject 
      DATA IQPML1/0/,IQQBL1/0/,IQQBCL1/0/,IQCDCL1/0/,IQPML0/0/
      DATA IQQBL0/0/,IQCDCL0/0/,IQPMH1/0/,IQQBH1/0/,IQQBCH1/0/
      DATA IQCDCH1/0/,IQPMH0/0/,IQQBH0/0/,IQCDCH0/0/
      data nreject/0/
      DATA FIRST/.TRUE./
      DATA NPRIN/5/
c      Data Fact/1.8d0/
      Data Fact/10d0/
      IDEBUGO = IDEBUG
   10 CONTINUE
      NIN=NIN+1
c      IF(MOD(NOUT,1).EQ.0) write(6,*) ' start: Nr of events ',NOUT
      IF(FIRST) THEN
         first = .FALSE.
         remPARJ32 = PARJ(32)
         IPROO = IPRO
         IDIRO = IDIR
         NGO = NG
         NPOMO = NPOM
         KPAO = KPA
         KINT22 = KINT(2,2)
         FIRST=.FALSE.
         MP = P(2,5)
         T2MAXO = T2MAX
         NDIS=0
         NQPMDS=0
         NQQBDS=0
         NQQBCDS=0
         NQQBBDS=0
         NQCDCDS=0
         NDIF=0
         NQPMDF=0
         NQQBDF=0
         NQQBCDF=0
         NQQBBDF=0
         NQCDCDF=0
         NPI=0
         NQPMPI=0
         NQQBPI=0
         NQQBCPI=0
         NQQBBPI=0
         NQCDCPI=0
         NRPA=0
         NRPB=0
         NRPC=0
         NRPD=0
         NRPE=0
         NRPF=0
         NRPG=0
         NPRT=0
         IHERPYSO = IHERPYS
         IHFO = IHF
      ELSE
         PARJ(32)=remPARJ32
         IHERPYS = IHERPYSO
         IHF = IHFO
         IPRO = IPROO
         IDIR=IDIRO
         NG = NGO
         NPOM = NPOMO
         KPA = KPAO
         KINT(2,2) = KINT22
      ENDIF
      T2GKI = 0.
      XFGKI = 0.
      PT2H = 0.
      SHH = 0.
      IWEI = 0
      IGENFL = 0
      IMIX = 0
      DO 20 I=1,20
         XV(I)= 0.D0
   20 X(I) = 0.D0
      NDIMEN = NDIM
      Igenfl = 0
      IF(IDISDIF.GE.1) IHERPYS = 0
      CALL RAPGEN(NDIM,XG)
      DO 30 I=1,NDIM
         X(I) = XG(I)
   30 XV(I) = X(I)
      IDIRO = IDIR
      IF(IDISDIF.GE.1) THEN
         IGENFL = 0
         XPRT = xpr
         Q2T = sngl(q2q)
         IDIRO = IDIR
         IF(XPRT.LT.XX(1).OR.XPRT.GT.XX(NBX) .OR.Q2T.LT.Q2X(1)
     +   .OR.Q2T.GT.Q2X(NBQ2)) THEN
            WRITE(6,*) ' EVENT: X or Q2 values outside grid '
            WRITE(6,*) ' X_min ',XX(1),' X_max ',XX(NBX),
     +         ' actual X ', XPRT
            WRITE(6,*) ' Q2_min ',Q2X(1),' Q2_max ',Q2X(NBQ2), ' '
     +         //'actual Q2 ',Q2T
            IF(XPRT.LT.XX(1)) XPRT=XX(1)
            IF(XPRT.GT.XX(NBX)) XPRT=XX(NBX)
            IF(Q2T.LT.Q2X(1)) Q2T = Q2X(1)
            IF(Q2T.GT.Q2X(NBQ2)) Q2T = Q2X(NBQ2)
         ENDIF
         IX = 0
  130    IX = IX + 1
         IF(XPRT.GT.XX(IX+1)) GOTO 130
         IQ = 0
  140    IQ = IQ + 1
         IF(Q2T.GT.Q2X(IQ+1)) GOTO 140

         CALL rgstfu(2212,XPRT,Q2T,XPQ2)
         IHERPYS = 1
         CALL rgstfu(2212,XPRT,Q2T,XPQ2_DIFF)
         idir = IDIRO
         IHERPYS = IHERPYSO
         nflavt = K(nia2,2)
         if(nflavt.eq.21) nflavt = 0
         
         F2DIS = xpq2(nflavt)
         F2DIF = xpq2_diff(nflavt)
         if(F2DIF.gt.F2DIS) F2DIS=F2DIF+0.001*F2DIF
         
         F2PI = 0
         RF2DIS = (F2DIS-F2DIF-F2PI)/F2DIS

         IF(RF2DIS.LE.0.0D0) THEN
            write(6,*) ' pomeron contribution larger than '
     +      //'total'
            write(6,*) ' TOT= ',F2DIS,' DIF= ',F2DIF,' PION= ',
     +      F2PI
            write(6,*) ' non-dffraction forced ..... mixing might be wrong '
            stop

         ENDIF
         RF2DIF = F2DIF/F2DIS
         RF2PI = F2PI/F2DIS
         KPA = KPAO
         RNTEST = draprn()
         IF(RNTEST.LT.RF2DIS) THEN
            IDIR = 1
		IDIRINT = IDIR
            KINT(2,2)=2212
            NDIMEN = NDIM
            WTMA = 999999.D0
	      NDIS = NDIS +1 
         ELSEIF(RNTEST.LT.(RF2DIS+RF2DIF)) THEN
            NTRY = 0
            NG = NGO
            NPOM =NPOMO
            IDIR = 0
            KINT(2,2) = 990
            WTMA = POMAX(IX,IQ)
            if(wtma.le.0) then
              write(6,*) ' eventpp: wtma < 0'
              goto 10
            endif  
            NDIMEN = NDIM + 2
c generate xr and t2
  150       CALL draprnV(XG,2)
            NTRY = NTRY + 1
            DO 160 I=1,2
               XV(NDIM+I)=XG(I)
  160       X(NDIM+I)=XG(I)
            XMAX = 1.D0 - XF
            XP2T = XPRT
            IF(XMAX.LE.XP2T) THEN
               write(6,*) '3. xmax<xp2t ',xmax,xp2t
               GOTO 10
            endif
            IF(NTRY.GT.NTRMAX) THEN
               if(nreject.le.1) then
               write(6,*) ' too many trials',xprt,t2,xr,wtdist,wtt,wpa
               write(6,*) ' too many trials for DIF;',
     +          ' event rejected,mix,po',ntry
               endif
               nreject = nreject + 1
               GOTO 10
            ENDIF
            XR = XP2T*(XMAX/XP2T)**X(NDIM+1)
            T2MIN = MP*MP*XR*XR/(1.D0-XR)
            T2MAX = T2MAXO

            T2 = T2MIN*((T2MAX/T2MIN)**X(NDIM+2))
            WTDIST = 0.D0
            CALL RAT2DI(KINT(2,2),XR,-T2,WTDIST)
            WTDIST = WTDIST * T2*DLOG(T2MAX/T2MIN)*XR*DLOG(XMAX/XP2T)
c include that when ng <0 or ng>30 because weight is only included in RASTFU,
c in RAT2DI weight=1.
            XPR = SNGL(XP2T)
            T2GKI = SNGL(-T2)
            XFGKI = SNGL(XR)
            CALL RASTFU(990,SNGL(XP2T/XR),Q2T,XPQ)
            if(NFLAVT.gt.NFLAV) then
               write(6,*) 'eventpp NFLAVT ',NFLAVT
            endif
            WTDIST = WTDIST*DBLE(XPQ(NFLAVT))

c add safety margin Fact for rejection
            WTT = WTDIST/(Fact*WTMA)
            if(WTT.GT.1) then
               write(6,*) ' eventpp: diffraction FATAL : WTT > 1',
     &           Wtt,Wtdist,WTMA
            endif
            WTRN = draprn()
            IF(WTRN.GT.WTT) THEN
               GOTO 150
            ENDIF
            WPA = FXN1(X,WDUM)
            IF(WPA.LE.0.0) then
                GOTO 150
            endif
            NDIF = NDIF + 1
         ELSEIF(RNTEST.LT.(RF2DIS+RF2DIF+RF2PI)) THEN
            NTRY = 0
            IDIR = 0
            NG = 20
            NPOM = 20
            KINT(2,2)=211
            WTMA = PIMAX(IX,IQ)
            NDIMEN = NDIM + 2
c generate xr and t2
  170       CALL draprnV(XG,2)
            DO 180 I=1,2
               XV(NDIM+I)=XG(I)
  180       X(NDIM+I)=XG(I)
            NTRY = NTRY + 1
            XMAX = 1.D0 - XF
            XP2T = XPRT
            IF(XMAX.LE.XP2T) THEN
               write(6,*) '4. xmax<xp2t ',xmax,xp2t
               GOTO 10
            endif
            IF(NTRY.GT.NTRMAX) THEN
               write(6,*) ' too many trials for PI;',
     +         '  event rejected,mix,pi'
               GOTO 10
            ENDIF
            XR = XP2T*(XMAX/XP2T)**X(NDIM+1)
            T2MIN = MP*MP*XR*XR/(1.D0-XR)
            T2MAX = T2MAXO

            T2 = T2MIN*((T2MAX/T2MIN)**X(NDIM+2))
            WTDIST = 0.D0
            CALL RAT2DI(KINT(2,2),XR,-T2,WTDIST)
            WTDIST = WTDIST * T2*DLOG(T2MAX/T2MIN)*XR*DLOG(XMAX/XP2T)
c include that when ng <0 or ng>30 because weight is only included in RASTFU,
c in RAT2DI weight=1.
            XPR = SNGL(XP2T)
            T2GKI = SNGL(-T2)
            XFGKI = SNGL(XR)
            CALL RASTFU(990,SNGL(XP2T/XR),Q2T,XPQ)
            WTDIST = WTDIST*DBLE(XPQ(KPA))
c add safety margin fact for rejection
            WTT = WTDIST/(Fact*WTMA)
            WTRN = draprn()
            IF(WTRN.GT.WTT) THEN
               GOTO 170
            ENDIF
            WPA = FXN1(X,WDUM)
            IF(WPA.LE.0.0) GOTO 170
            NPI = NPI  + 1
         ELSE
            WRITE(6,*) ' something"s wrong here '
            WRITE(6,*) ' RNTEST = ',RNTEST,' RF2DIS = ',RF2DIS, ' '
     +      //'RF2DIF = ',RF2DIF,' RF2PI = ',RF2PI
         ENDIF
         WPA = FXN1(X,WDUM)
      ENDIF
      IHERPYS =0

      DO 190 I=1,NDIMEN
  190 XV(I) = X(I)
      IGENFL = 1
      WPA = FXN1(X,WDUM)
      IF(WPA.LE.0.0) THEN
         GOTO 10
      ENDIF
      IDEBUG = IDEBUGO
      IGENFL = 1
c zero MEINFO
      ZQGKI = 99999.
      XPGKI = 99999.
      PHITGKI = 99999.
c............
  200 CONTINUE
      IWEI = 0


      IF(IFPS.NE.10) THEN
         IHF=0
         IF(IPRO.EQ.99) THEN
         ELSE
            CALL PPMEPS
         ENDIF
         IHF = IHFO
c check on ppmeps and pyremn
         IF(LST(21).NE.0) THEN
            LERR(LST(21)) = LERR(LST(21)) + 1
            GOTO 10
         ENDIF
      ELSE
         CALL PPMEPS
c check on ppmeps and pyremn
         IF(LST(21).NE.0) THEN
            LERR(LST(21)) = LERR(LST(21)) + 1
            GOTO 10
         ENDIF
         MSTA(5)=0
         CALL AREXEC
      ENDIF
      IF(NFRAG.EQ.10.AND.NPOM.NE.20.AND.NPOM.NE.21
     +           .AND.IDIR.EQ.0) CALL PRODIFF

      IF(NFRAG.GE.1.and.ILHA.eq.0) CALL pyexec
      CALL pyrobo(0,N,0.d0,0.d0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
c only for decay of  p - diss
      IF(NFRAG.EQ.10.AND.NPOM.NE.20.AND.NPOM.NE.21
     +           .AND.IDIR.EQ.0) CALL PDISDC
c....
      IF(MSTU(24).NE.0) THEN
         WRITE(6,*) 'MSTU(24)= ',MSTU(24)
         CALL pylist(1)
      ENDIF
      NOUT = NOUT + 1
      IF(MOD(NOUT,1000).EQ.0) then 
         write(6,*) ' Nr of events ',NOUT,' non-diff ',NDIS,' diff
     &    ',NDIF,nreject
         
      endif 
      IF(IPRO.EQ.18) THEN
         NPRT = NPRT + 1
         IF(IRESPRO.EQ.1) NRPA=NRPA+1
         IF(IRESPRO.EQ.2) NRPB=NRPB+1
         IF(IRESPRO.EQ.3) NRPC=NRPC+1
         IF(IRESPRO.EQ.4) NRPD=NRPD+1
         IF(IRESPRO.EQ.5) NRPE=NRPE+1
         IF(IRESPRO.EQ.6) NRPF=NRPF+1
         IF(IRESPRO.EQ.7) NRPG=NRPG+1
         IF(IRESPRO.EQ.8) NRPH=NRPH+1
      ENDIF
      IF(NOUT.LE.NPRIN) THEN
         write(6,*) ' event: IPRO = ',IPRO
         CALL pylist(1)
      ENDIF
      NG = NGO
      NPOM = NPOMO
      RETURN
      END
