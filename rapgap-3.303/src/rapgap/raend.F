      SUBROUTINE RAEND(IK)
      IMPLICIT None
      CHARACTER* 16 CNAM
      CHARACTER* 6 CBEAM1,CBEAM2
      CHARACTER * 7 CINT
	Integer NQPM,NQQB,NQQBC,NQQBB,NQCDC,IMIX,IK,LERR
	Double Precision WMAX
      COMMON /QCDEV/ NQPM,NQQB,NQQBC,NQQBB,NQCDC
      COMMON /OALPHAS/ WMAX,IMIX
      COMMON/ERR/LERR(100)
#include "rghs45.inc"
#include "rgnlo.inc"
#include "rgcolcon.inc"
      Integer NTHRE,NVMRE,IVM
      COMMON/HERTHET/ NTHRE,NVMRE
      COMMON/VMESON/IVM
	Integer NDIS,NQPMDS,NQQBDS,NQQBCDS,NQQBBDS,NQCDCDS,
     +NDIF,NQPMDF,NQQBDF,NQQBCDF,NQQBBDF,NQCDCDF,
     +NPI,NQPMPI,NQQBPI,NQQBCPI,NQQBBPI,NQCDCPI,
     +NRPA,NRPB,NRPC,NRPD,NRPE,NRPF,NRPG,NRPH,NPRT

      COMMON/NEVOUT/NDIS,NQPMDS,NQQBDS,NQQBCDS,NQQBBDS,NQCDCDS,
     +              NDIF,NQPMDF,NQQBDF,NQQBCDF,NQQBBDF,NQCDCDF,
     +              NPI,NQPMPI,NQQBPI,NQQBCPI,NQQBBPI,NQCDCPI,
     +              NRPA,NRPB,NRPC,NRPD,NRPE,NRPF,NRPG,NRPH,NPRT
      REAL ACC1,ACC2
	Integer IINT,NCB
      COMMON /INTEGR/ ACC1,ACC2,IINT,NCB
#include "rgpart.inc"
#include "rgeffic.inc"
#include "rgparas.inc"
#include "rglujets.inc"
#include "rgdisdif.inc"
#include "rgluco.inc"
      Double Precision SIGQPM,ESIGQPM,SIGQQB,ESIGQQB,SIGQQBC,ESIGQQBC
	Double Precision SIGQQBB,ESIGQQBB,SIGQCDC,ESIGQCDC,XVIS,XVMVIS
      IF(IK.EQ.10) WRITE(6,*) 'TIME LIMIT REACHED........'
      CALL PYNAME(KBEAM(1,2),CNAM)
      CBEAM1 = CNAM
      CALL PYNAME(KBEAM(2,2),CNAM)
      CBEAM2 = CNAM
      CALL PYNAME(KINT(2,2),CNAM)
      CINT = CNAM
      IF(IMIX.NE.1) THEN
         IF(IPRO.EQ.10.OR.IPRO.EQ.13) THEN
            WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> q q_'
     +      //'bar X"'
         ELSEIF(IPRO.EQ.11.OR.IPRO.EQ.14) THEN
            WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> Q Q_'
     +      //'bar X"'
         ELSEIF(IPRO.EQ.12) THEN
            WRITE(6,*) ' x - section for ',CBEAM1,CBEAM2,' --> q q_bar '
     +      //'X"'
         ENDIF
         WRITE(6,10000) AVGI,SD
      ENDIF
10000 FORMAT('  sigma = ',G16.7,' nb   +/- ',G16.7)
      IF(IMIX.NE.1) THEN
         IF(IPRO.EQ.10.OR.IPRO.EQ.13)
     +   WRITE(6,*) ' [ gamma g_',CINT,' --> q q_bar] ',avgi
         IF((IPRO.EQ.11.OR.IPRO.EQ.14).and.IHFLA.EQ.4)
     +   WRITE(6,*) ' [ gamma g_',CINT,' --> c c_bar ] ',avgi
         IF((IPRO.EQ.11.OR.IPRO.EQ.14).and.IHFLA.EQ.5)
     +   WRITE(6,*) ' [ gamma g_',CINT,' --> b b_bar ] ',avgi
         IF(IPRO.EQ.12) THEN
            WRITE(6,*) ' [ gamma q_',CINT,' --> q ] ',avgi
         ENDIF
         IF(IPRO.EQ.18) THEN
            IF(IK.EQ.20) THEN
               WRITE(6,*) ' QCD processes: '
               IF(IRPA.EQ.1) WRITE(6,1001)
     +         avgi*DFLOAT(NRPA)/DFLOAT(NPRT)
1001  Format(' [ g + g     --> q + q_bar ]...',E10.3)   
               IF(IRPB.EQ.1) WRITE(6,1002) ,
     +         avgi*DFLOAT(NRPB)/DFLOAT(NPRT)
1002  Format(' [ g + g     --> g + g]........',E10.3)    
               IF(IRPC.EQ.1) WRITE(6,1003) 
     +         avgi*DFLOAT(NRPC)/DFLOAT(NPRT)
1003  Format(' [ g + q     --> g + q ].......',E10.3)    
               IF(IRPD.EQ.1) WRITE(6,1004) 
     +         avgi*DFLOAT(NRPD)/DFLOAT(NPRT)
1004  Format(' [ q + q_bar --> g + g ].......',E10.3)    
               IF(IRPE.EQ.1) WRITE(6,1005) 
     +         avgi*DFLOAT(NRPE)/DFLOAT(NPRT)
1005  Format(' [ q + q_bar --> q + q_bar ]...',E10.3)    
               IF(IRPF.EQ.1) WRITE(6,1006) 
     +         avgi*DFLOAT(NRPF)/DFLOAT(NPRT)
1006  Format(' [ q + q     --> q + q ].......',E10.3)    
               IF(IRPG.EQ.1)
     +	   WRITE(6,1007) 
     +         avgi*DFLOAT(NRPG)/DFLOAT(NPRT)
1007  Format(' Mueller/Tang  [ q + q     --> q + q ] ',E10.3)    
               IF(IRPG.EQ.2)
     +	   WRITE(6,1008) 
     +         avgi*DFLOAT(NRPG)/DFLOAT(NPRT)
1008  Format(' massive glu   [ q + q     --> q + q ] ',E10.3)    
               IF(IRPH.EQ.1)
     +	   WRITE(6,1009) 
     +         avgi*DFLOAT(NRPH)/DFLOAT(NPRT)
1009  Format('  [ q + g     --> gamma + g ] ',E10.3)    
            ENDIF
         ENDIF
      ENDIF
      IF(IMIX.EQ.1) THEN
         IF(IDISDIF.EQ.0) THEN
            SIGQPM = AVGI*DFLOAT(NQPM)/DFLOAT(NOUT)
            ESIGQPM = AVGI*DSQRT(DFLOAT(NQPM))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma q_',CINT,' --> q ] ',SIGQPM,' +/- ',
     +      ESIGQPM
            SIGQQB = AVGI*DFLOAT(NQQB)/DFLOAT(NOUT)
            ESIGQQB= AVGI*DSQRT(DFLOAT(NQQB))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> q q_bar] ',SIGQQB,
     +      ' +/- ',ESIGQQB
            SIGQQBC = AVGI*DFLOAT(NQQBC)/DFLOAT(NOUT)
            ESIGQQBC = AVGI*DSQRT(DFLOAT(NQQBC))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> c c_bar] ',SIGQQBC,
     +      ' +/- ',ESIGQQBC
            SIGQQBB = AVGI*DFLOAT(NQQBB)/DFLOAT(NOUT)
            ESIGQQBB = AVGI*DSQRT(DFLOAT(NQQBB))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> b b_bar] ',SIGQQBB,
     +      ' +/- ',ESIGQQBB
            SIGQCDC = AVGI*DFLOAT(NQCDC)/DFLOAT(NOUT)
            ESIGQCDC = AVGI*DSQRT(DFLOAT(NQCDC))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma q_',CINT,' --> q g] ',SIGQCDC, ' +/- '
     +      ,ESIGQCDC
         ELSE
            WRITE(6,*) ' global x sections '
            SIGQPM = AVGI*DFLOAT(NQPM)/DFLOAT(NOUT)
            ESIGQPM = AVGI*DSQRT(DFLOAT(NQPM))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma q_',CINT,' --> q ] ',SIGQPM,' +/- ',
     +      ESIGQPM
            SIGQQB = AVGI*DFLOAT(NQQB)/DFLOAT(NOUT)
            ESIGQQB= AVGI*DSQRT(DFLOAT(NQQB))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> q q_bar] ',SIGQQB,
     +      ' +/- ',ESIGQQB
            SIGQQBC = AVGI*DFLOAT(NQQBC)/DFLOAT(NOUT)
            ESIGQQBC = AVGI*DSQRT(DFLOAT(NQQBC))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> c c_bar] ',SIGQQBC,
     +      ' +/- ',ESIGQQBC
            SIGQQBB = AVGI*DFLOAT(NQQBB)/DFLOAT(NOUT)
            ESIGQQBB = AVGI*DSQRT(DFLOAT(NQQBB))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> b b_bar] ',SIGQQBB,
     +      ' +/- ',ESIGQQBB
            SIGQCDC = AVGI*DFLOAT(NQCDC)/DFLOAT(NOUT)
            ESIGQCDC = AVGI*DSQRT(DFLOAT(NQCDC))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma q_',CINT,' --> q g] ',SIGQCDC, ' +/- '
     +      ,ESIGQCDC
            WRITE(6,*) ' DIS x section: ',AVGI*DFLOAT(NDIS)/
     +      DFLOAT(NOUT)
            SIGQPM = AVGI*DFLOAT(NQPMDS)/DFLOAT(NOUT)
            ESIGQPM = AVGI*DSQRT(DFLOAT(NQPMDS))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma q_',CINT,' --> q ] ',SIGQPM,' +/- ',
     +      ESIGQPM
            SIGQQB = AVGI*DFLOAT(NQQBDS)/DFLOAT(NOUT)
            ESIGQQB= AVGI*DSQRT(DFLOAT(NQQBDS))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> q q_bar] ',SIGQQB,
     +      ' +/- ',ESIGQQB
            SIGQQBC = AVGI*DFLOAT(NQQBCDS)/DFLOAT(NOUT)
            ESIGQQBC = AVGI*DSQRT(DFLOAT(NQQBCDS))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> c c_bar] ',SIGQQBC,
     +      ' +/- ',ESIGQQBC
            SIGQQBB = AVGI*DFLOAT(NQQBBDS)/DFLOAT(NOUT)
            ESIGQQBB = AVGI*DSQRT(DFLOAT(NQQBBDS))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> b b_bar] ',SIGQQBB,
     +      ' +/- ',ESIGQQBB
            SIGQCDC = AVGI*DFLOAT(NQCDCDS)/DFLOAT(NOUT)
            ESIGQCDC = AVGI*DSQRT(DFLOAT(NQCDCDS))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma q_',CINT,' --> q g] ',SIGQCDC, ' +/- '
     +      ,ESIGQCDC
            WRITE(6,*) ' DIF x section ',AVGI*DFLOAT(NDIF)/DFLOAT(NOUT)
            SIGQPM = AVGI*DFLOAT(NQPMDF)/DFLOAT(NOUT)
            ESIGQPM = AVGI*DSQRT(DFLOAT(NQPMDF))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma q_',CINT,' --> q ] ',SIGQPM,' +/- ',
     +      ESIGQPM
            SIGQQB = AVGI*DFLOAT(NQQBDF)/DFLOAT(NOUT)
            ESIGQQB= AVGI*DSQRT(DFLOAT(NQQBDF))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> q q_bar] ',SIGQQB,
     +      ' +/- ',ESIGQQB
            SIGQQBC = AVGI*DFLOAT(NQQBCDF)/DFLOAT(NOUT)
            ESIGQQBC = AVGI*DSQRT(DFLOAT(NQQBCDF))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> c c_bar] ',SIGQQBC,
     +      ' +/- ',ESIGQQBC
            SIGQQBB = AVGI*DFLOAT(NQQBBDF)/DFLOAT(NOUT)
            ESIGQQBB = AVGI*DSQRT(DFLOAT(NQQBBDF))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma g_',CINT,' --> b b_bar] ',SIGQQBB,
     +      ' +/- ',ESIGQQBB
            SIGQCDC = AVGI*DFLOAT(NQCDCDF)/DFLOAT(NOUT)
            ESIGQCDC = AVGI*DSQRT(DFLOAT(NQCDCDF))/DFLOAT(NOUT)
            WRITE(6,*) ' [ gamma q_',CINT,' --> q g] ',SIGQCDC, ' +/- '
     +      ,ESIGQCDC
            IF(IDISDIF.EQ.2) THEN
               WRITE(6,*) ' PI x section ',AVGI*DFLOAT(NPI)/
     +         DFLOAT(NOUT)
               SIGQPM = AVGI*DFLOAT(NQPMPI)/DFLOAT(NOUT)
               ESIGQPM = AVGI*DSQRT(DFLOAT(NQPMPI))/DFLOAT(NOUT)
               WRITE(6,*) ' [ gamma q_',CINT,' --> q ] ',SIGQPM,' +/- '
     +         ,ESIGQPM
               SIGQQB = AVGI*DFLOAT(NQQBPI)/DFLOAT(NOUT)
               ESIGQQB= AVGI*DSQRT(DFLOAT(NQQBPI))/DFLOAT(NOUT)
               WRITE(6,*) ' [ gamma g_',CINT,' --> q q_bar] ',SIGQQB,
     +         ' +/- ',ESIGQQB
               SIGQQBC = AVGI*DFLOAT(NQQBCPI)/DFLOAT(NOUT)
               ESIGQQBC = AVGI*DSQRT(DFLOAT(NQQBCPI))/DFLOAT(NOUT)
               WRITE(6,*) ' [ gamma g_',CINT,' --> c c_bar] ',SIGQQBC,
     +         ' +/- ',ESIGQQBC
               SIGQQBB = AVGI*DFLOAT(NQQBBPI)/DFLOAT(NOUT)
               ESIGQQBB = AVGI*DSQRT(DFLOAT(NQQBBPI))/DFLOAT(NOUT)
               WRITE(6,*) ' [ gamma g_',CINT,' --> b b_bar] ',SIGQQBB,
     +         ' +/- ',ESIGQQBB
               SIGQCDC = AVGI*DFLOAT(NQCDCPI)/DFLOAT(NOUT)
               ESIGQCDC = AVGI*DSQRT(DFLOAT(NQCDCPI))/DFLOAT(NOUT)
               WRITE(6,*) ' [ gamma q_',CINT,' --> q g] ',SIGQCDC,
     +         ' +/- ',ESIGQCDC
            ENDIF
         ENDIF
      ENDIF
      IF(IK.EQ.20) THEN
         IF(IINT.EQ.0.and.IHERAC.EQ.0)  CALL SPINFO( 6 )
         IF(IHERAC.EQ.1) THEN
            write(6,*) ' Nr of events rej.by theta cut ',NTHRE
            XVIS = AVGI*DFLOAT(NOUT)/DFLOAT(NOUT+NTHRE)
            write(6,*) ' visible x section = ',XVIS
            IF(IVM.NE.0) THEN
               XVMVIS = AVGI*DFLOAT(NOUT)/DFLOAT(NIN)
               write(6,*) ' vector meson x section = ',XVMVIS
            ENDIF
            AVGI = XVIS
         ENDIF
         WRITE(6,*) ' Nr of events generated : ',NIN
         WRITE(6,*) ' Nr of events written: ',NOUT
         WRITE(6,*) ' Error summary on event generation '
         WRITE(6,*) ' Errors and their meaning:'
	   if(IBS.EQ.1) then
            WRITE(6,*) ' ELQQFNL: neagtive ME......: ',NEGNLO	   
            WRITE(6,*) ' ELQQFNL: out of NTOTNLO...: ',NTOTNLO	   
	   endif
         WRITE(6,*) ' QCDMIX: WPA>WMAX..........: ',LERR(30)
         WRITE(6,*) ' QCDMIX: WPA<0.0...........: ',LERR(31)
         WRITE(6,*) ' QCDMIX: NCALL>1000........: ',LERR(32)
         WRITE(6,*) ' meps x > 0.999............: ',LERR(45)
         WRITE(6,*) ' meos boost PS error.......: ',LERR(46)
         WRITE(6,*) ' meps energy not conserved.: ',LERR(100)
         WRITE(6,*) ' rgremn frag. cuts.........: ',LERR(48),LERR(49)
         WRITE(6,*) ' PYPREP error..............: ',LERR(50)
         WRITE(6,*) ' rasspa check 1st..........: ',LERR(55)
         WRITE(6,*) ' rasspa no more memory.....: ',LERR(51),LERR(52)
         WRITE(6,*) ' rasspa boost error........: ',LERR(53)
         WRITE(6,*) ' rasspa xfb(iflb)=0........: ',LERR(54)
         WRITE(6,*) ' rasspa NTRY > 5000........: ',LERR(56)
         WRITE(6,*) ' rasspa XB<0...............: ',LERR(57)

         IF(IPRO.GT.1000) CALL HERACL(AVGI,SD,3)

      ENDIF
      RETURN
      END
