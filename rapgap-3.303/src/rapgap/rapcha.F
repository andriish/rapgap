
      SUBROUTINE RAPCHA

*#**********************************************************************
*#
*#    SUBROUTINE RAPCHA
*#
*# PURPOSE: TO CHANGE THE DEFAULT
*#          PARAMETERS OF RAPGAP.
*#
*#**********************************************************************
	Integer nloop,nmax
	Parameter (nmax=1000)
      character *4   para
      character*50 Cdef,Cea
	Real Rval 
	Integer Ival,Ld,Le
      character*132 Cval 
	Common/steering/Nloop,Ld(nmax),Le(nmax),Ival(nmax),
     &	para(nmax),Rval(nmax),Cval(nmax)

      CHARACTER*4     PARAM
*
#include "rgpypars.inc"
#include "rgpara.inc"
#include "rgluco.inc"
#include "rgpart.inc"
#include "rgdiffr.inc"
#include "rgdisdif.inc"
#include "rgfull.inc"
#include "rgscq2.inc"
#include "rglq2.inc"
#include "pshwr.inc"
#include "rgprkt.inc"
#include "rgpqcdpom.inc"
#include "rgnlo.inc"
#include "rgpdiss.inc"
#include "rgherup.inc"
#include "rgcolcon.inc"
      character*512 PDFPATH
      common/RAPPDF/PDFPATH

      INTEGER IGROUP1, ISET1,idif1
	Common/strnew/IGROUP1,ISET1,idif1
      COMMON/VMESON/IVM
      DOUBLE PRECISION C1,Cg
      COMMON/BUCHMUE/C1,Cg
      REAL ACC1,ACC2
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
      Integer IGLU
      Common/gluon/Iglu
      INTEGER LST,IRES
      COMMON/EPPARA/LST(30),IRES(2)

C---Determine the machine accuacy EPSMAC (D= 2.E-06)
      EPSMAC = 0.5
      DO 10 I= 1, 100
         EPSMAC = EPSMAC * 0.5
         IF ((1.0+EPSMAC) .EQ. 1.0) GO TO 20
   10 CONTINUE
      EPSMAC = 1.0E-6
   20 EPSMAC = 2.0 * EPSMAC
      WRITE(6,10000) EPSMAC
10000 FORMAT(1X,/,' COMPUTED MACHINE ACCURACY = ',G10.2,/)


C---Interpret the text containing the parameters to be changed
      If(Nloop.ge.1) Then
            
		Do 30 I=1,Nloop

            PARAM = Para(I)
            INT = Ival(I)
            REA = Rval(I)
            CEA = Cval(I)

c		write(6,*) ' params', param,int,rea,cea
            IF (PARAM.EQ.'LEPI') THEN
               IDEF = KE
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  KE = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'HADI') THEN
               IDEF = KP
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  KP = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'KBE1') THEN
               IDEF = KE
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  KE = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRE1') THEN
               IDEF = Ires(1)
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ires(1) = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'KBE2') THEN
               IDEF = KP
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  KP = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRE2') THEN
               IDEF = Ires(2)
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ires(2) = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'NFRA') THEN
               IDEF = NFRAG
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NFRAG = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IFPS') THEN
               IDEF = IFPS
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IFPS = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
             IF (PARAM.EQ.'IPST') THEN
               IDEF = IPST
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IPST = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
           IF (PARAM.EQ.'IRAS') THEN
               IDEF = IRUNA
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRUNA = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRAM') THEN
               IDEF = IRUNAEM
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRUNAEM = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IQ2S') THEN
               IDEF = IQ2
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IQ2 = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IPRO') THEN
               IDEF = IPRO
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IPRO = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'INGG') THEN
               IDEF = NG
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NG = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'INPR') THEN
               IDEF = MSTP(51)
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  MSTP(51) = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'INPI') THEN
               IDEF = MSTP(52)
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  MSTP(52) = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'INGA') THEN
               IDEF = MSTP(55)
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  MSTP(55) = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IDIR') THEN
               IDEF = IDIR
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IDIR = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IAKT') THEN
               IDEF = IALMKT
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IALMKT = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IDIF') THEN
               IDEF = IDISDIF
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IDISDIF = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IFUL') THEN
               IDEF = IFULL
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IFULL = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IGRI') THEN
               IDEF = IQCDGRID
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IQCDGRID = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IPOM') THEN
               IDEF = NPOM
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NPOM = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IVME') THEN
               IDEF = IVM
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IVM = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IREM') THEN
               IDEF = IREM
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IREM = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'NFLA') THEN
               IDEF = NFLAV
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NFLAV = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'NFQC') THEN
               IDEF = NFLQCDC
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NFLQCDC = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'INTE') THEN
               IDEF = INTER
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  INTER = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'IHFL') THEN
               IDEF = IHFLA
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IHFLA = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'IRPA') THEN
               IDEF = IRPA
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPA = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPB') THEN
               IDEF = IRPB
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPB = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPC') THEN
               IDEF = IRPC
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPC = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPD') THEN
               IDEF = IRPD
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPD = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPE') THEN
               IDEF = IRPE
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPE = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPF') THEN
               IDEF = IRPF
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPF = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPG') THEN
               IDEF = IRPG
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPG = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPH') THEN
               IDEF = IRPH
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPH = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRPI') THEN
               IDEF = IRPI
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRPI = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'ISET') THEN
               IDEF = ISET
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ISET = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IP2 ') THEN
               IDEF = IP2
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IP2 = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF

            IF (PARAM.EQ.'ISEH') THEN
               IDEF = ISEMIH
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ISEMIH = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IORD') THEN
               IDEF = IORDER
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IORDER = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ICCF') THEN
               IDEF = ICCFM
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ICCFM = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IALP') THEN
               IDEF = IALPS
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IALPS = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ITIM') THEN
               IDEF = ITIMSHR
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ITIMSHR = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ISOG') THEN
               IDEF = ISOFTG
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ISOFTG = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'NCAL') THEN
               IDEF= NCB
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  NCB = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRSP') THEN
               IDEF= IRspl
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IRspl = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'Iqqg') THEN
               IDEF= Iqqg
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Iqqg = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'Igkt') THEN
               IDEF= Igamkt
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Igamkt = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IGLU') THEN
               IDEF = IGLU
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IGLU = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IBS ') THEN
               IDEF = IBS
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IBS = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IJET') THEN
               IDEF = ijetset
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ijetset = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IFRP') THEN
               IDEF = ifragp
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ifragp = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IRGL') THEN
               IDEF = i_glu
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  i_glu = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ILHA') THEN
               IDEF = Ilha
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  Ilha = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'IGR ') THEN
               IDEF = IGROUP1
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  IGROUP1 = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ISE ') THEN
               IDEF = ISET1
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ISET1 = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
             IF (PARAM.EQ.'ILAM') THEN
               IDEF = ILAMQCD
               IF (ABS(IDEF-INT).GT.EPSMAC) THEN
                  ILAMQCD = INT
                  PRINT 10100, PARAM,LND,IDEF,INT
               ELSE
                  PRINT 10200, PARAM,LND,IDEF,INT
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PBE2') THEN
               RDEF = PPIN
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PPIN = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PBE1') THEN
               RDEF = PLEPIN
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PLEPIN = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF

            IF (PARAM.EQ.'PPIN') THEN
               RDEF = PPIN
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PPIN = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PLEP') THEN
               RDEF = PLEPIN
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PLEPIN = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'RRN2') THEN
               RDEF = SNGL(RN2)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  RN2 = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'EESP') THEN
               RDEF = SNGL(EPSP)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  EPSP = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'QMIN') THEN
               RDEF = SNGL(QMI)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  QMI = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'QMAX') THEN
               RDEF = SNGL(QMA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  QMA = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'YMIN') THEN
               RDEF = SNGL(YMI)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  YMI = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'YMAX') THEN
               RDEF = SNGL(YMA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  YMA = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'RALP') THEN
               RDEF = SNGL(ALPHP)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  ALPHP = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'RRXF') THEN
               RDEF = SNGL(XF)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  XF = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'RT2M') THEN
               RDEF = SNGL(T2MAX)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  T2MAX = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PTCU') THEN
               IPROF = IPRO
               IF(IPRO.EQ.12.OR.IPRO.EQ.1200) THEN
                  PT2CUT(10) = DBLE(REA)
                  PT2CUT(13) = DBLE(REA)
                  PT2CUT(15) = DBLE(REA)
                  PT2CUT(18) = DBLE(REA)
               ELSE
                  IF(IPRO.GT.1000) IPROF = IPRO/1000
                  IF(IPRO.EQ.14.OR.IPRO.EQ.11) THEN
                     PT2CUT(IPROF) = 0.D0
                  ELSE
                     RDEF = SNGL(PT2CUT(IPROF))
                     IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                        PT2CUT(IPROF) = DBLE(REA)
                        PRINT 10300, PARAM,LND,RDEF,REA
                     ELSE
                        PRINT 10400, PARAM,LND,RDEF,REA
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            IF (PARAM.EQ.'THEM') THEN
               RDEF = SNGL(THEMA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  THEMA = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'THEI') THEN
               RDEF = SNGL(THEMI)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  THEMI = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF

            IF (PARAM.EQ.'WEQ2') THEN
               Q2START = 0.D0
               Q2 = DBLE(REA)
               IF(Q2.GT.0.1D0) THEN
                  Q2START=Q2
                  WRITE(6,*)' GRACHA: Weighting below Q2=',Q2START
               ENDIF
            ENDIF

            IF (PARAM.EQ.'OMEG') THEN
               RDEF = SNGL(OMEG2)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  OMEG2 = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF

            IF (PARAM.EQ.'SCQ2') THEN
               RDEF = SNGL(SCALQ2)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  SCALQ2 = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF

            IF (PARAM.EQ.'LQ2S') THEN
               RDEF = SNGL(Q2SUPP)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  Q2SUPP = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'KT1 ') THEN
               RDEF = SNGL(PRKT1)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKT1 = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'KTM ') THEN
               RDEF = SNGL(PKTM)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PKTM = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'KT2 ') THEN
               RDEF = SNGL(PRKT2)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  PRKT2 = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'SCAL') THEN
               RDEF = SNGL(SCALFA)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  SCALFA = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'C1  ') THEN
               RDEF = SNGL(C1)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  C1 = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'Cg  ') THEN
               RDEF = SNGL(Cg)
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  Cg = DBLE(REA)
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ACC1') THEN
               RDEF = ACC1
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  ACC1 = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'ACC2') THEN
               RDEF = ACC2
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  ACC2 = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PEPS') THEN
               RDEF = peps
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  peps = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PRN2') THEN
               RDEF = prn2
               IF (ABS(RDEF-REA).GT.EPSMAC) THEN
                  prn2 = REA
                  PRINT 10300, PARAM,LND,RDEF,REA
               ELSE
                  PRINT 10400, PARAM,LND,RDEF,REA
               ENDIF
            ENDIF
            IF (PARAM.EQ.'PDF1') THEN
               CDEF = pdfname
c               write(6,*) ' PDFname1 ',pdfname,cea,mstp(51)
               IF (CDEF.NE.CEA) THEN
                  pdfname = CEA
                  PRINT 10500, PARAM,LND,CDEF,CEA
               ELSE
                  PRINT 10600, PARAM,LND,CDEF,CEA
               ENDIF
               if(pdfname.ne.' ') mstp(51)=999999
            ENDIF
            IF (PARAM.EQ.'PDF2') THEN 
               CDEF = phname
c               write(6,*) ' PDFname1 ',pdfname,cea,mstp(51)
               IF (CDEF.NE.CEA) THEN
                  phname = CEA
                  PRINT 10500, PARAM,LND,CDEF,CEA
               ELSE
                  PRINT 10600, PARAM,LND,CDEF,CEA
               ENDIF
               if(phname.ne.' ') mstp(53)=999999
            ENDIF
            IF (PARAM.EQ.'PDF3') THEN
               CDEF = piname
c               write(6,*) ' PDFname1 ',pdfname,cea,mstp(51)
               IF (CDEF.NE.CEA) THEN
                  piname = CEA
                  PRINT 10500, PARAM,LND,CDEF,CEA
               ELSE
                  PRINT 10600, PARAM,LND,CDEF,CEA
               ENDIF
               if(piname.ne.' ') mstp(52)=999999
            ENDIF
            IF (PARAM.EQ.'UPDF') THEN
               CDEF = PDFPATH
               write(6,*) ' test ',cea
               IF (CDEF.NE.CEA) THEN
                  PDFPATH = CEA
                  PRINT 10500, PARAM,LND,CDEF,CEA
               ELSE
                  PRINT 10600, PARAM,LND,CDEF,CEA
               ENDIF
            ENDIF

   30    CONTINUE
      ENDIF

10100 FORMAT(' rapcha: ',A4,'(',I4,')',
     +       ' changed from the default ',I8,' TO',I8)
10200 FORMAT(' rapcha (WARNING): ',A4,'(',I4,')',
     +       ' not changed from the default',I8,' TO',I8)
10300 FORMAT(' rapcha:',A4,'(',I4,')',
     +       ' changed from the default ',E10.3,' TO',E10.3)
10400 FORMAT(' rapcha (WARNING):',A4,'(',I4,')',
     +       ' not changed from the default ',E10.3,' TO',E10.3)

10500 FORMAT(' rapcha: ',A4,'(',I4,')',
     +       ' changed from the default ',A50,' TO ',A50)
10600 FORMAT(' rapcha (WARNING): ',A4,'(',I4,')',
     +       ' not changed from the default ',A50,' TO ',A50)


      RETURN
      END
