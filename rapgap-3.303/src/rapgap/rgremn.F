C **********************************************************************
      SUBROUTINE rgremn(IPU1,IPU2)
      IMPLICIT None
C...ADDS ON TARGET REMNANTS (ONE OR TWO FROM EACH SIDE) AND
C...INCLUDES PRIMORDIAL KT.
      DOUBLE PRECISION DBETAX,DBETAZ,DROBO(5)
      REAL CUT,XLP,YLP,W2LP,Q2LP,ULP
      DOUBLE PRECISION PARL
      Integer LLST
      COMMON /RAPTOU/ CUT(14),LLST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
#include "rglujets.inc"
#include "rgludat1.inc"
#include "rgherup.inc"
#include "rgrapgki.inc"
#include "rgluco.inc"

#include "rgparas.inc"
#include "rgpara1.inc"
c.hju
#include "rgprkt.inc"
#include "rgdisdif.inc"
      Double Precision PYPAR,PYVAR
      Integer IPY
      COMMON /RGPYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      REAL X,SH,TH,UH,Q2PY
      DOUBLE PRECISION ROBO,SINTH,PHIPT
      Integer ISUB,KFLCH,KFLSP,IS,ISN,KFL
      COMMON /PYPROC/ ISUB,KFL(3,2),X(2),SH,TH,UH,Q2PY
      Double Precision PMS,PSYS
      DIMENSION KFLCH(2),KFLSP(2),PMS(0:6),IS(2),ROBO(5)

c.hju
      DIMENSION PSYS(0:2,5),ISN(2)

      Integer LST,IRES
      COMMON/EPPARA/LST(30),IRES(2)
      REAL GGPT,GPPT
      COMMON /GRAKT/ GGPT,GPPT
      real prkt
      Double Precision pktmin,pktmax
      REAL PTSPL,PHISPL,CHI(2)
      Double Precision zz,pmst
      REAL SNGL
      Integer I,J,JT,pycomp,IFLS,IPU,IQ,IPU1,IPU2,ILEP,IP,NS,NTRY
      Integer IL,IR,IG,IMAX,IMIN
      Double Precision PEI,PE,PZI,PZ,SHS,SHR,PMMIN0,PMMIN1
      Double Precision PMMIN2,PMMIN,PPB,PNB,PMTB,PMTR,PMTL,PW1
      Double Precision SQLAM,SQSGN,RKR,RKL,BER,BEL,PZF,PT2,RQP
      Double Precision PEM,PZM,PDEV
      EXTERNAL pycomp
      Double Precision draprn
	Integer icallz
      External draprn
      DATA IPU,IQ/0,0/,PEI,PE,PZI,PZ,SHS/5*0./
c... hju data added
c      write(6,*) 'IRES(2) ',IRES(2)
c...      IRES(2) = 1
c... hju
C...FIND EVENT TYPE, SET POINTERS
c      write(6,*) 'rgremn test; T2gki,xfgki ',T2GKI,XFGKI
c      write(6,*) 'rgremn test'
      GGPT = -9999.
      GPPT = -9999.
      icallz = 0
      IF(IPU1.EQ.0.AND.IPU2.EQ.0) RETURN
      ILEP=0
      IF(IPU1.EQ.0) ILEP=1
      IF(IPU2.EQ.0) ILEP=2
      IF(ISUB.EQ.7) ILEP=-1
      IF(ILEP.EQ.1) IQ=21
      IF(ILEP.EQ.2) IQ=23
      IF(ILEP.eq.1.and.ILEPTO.eq.0) ILEP=0

      IF(ILEP.EQ.1.OR.(IRES(1).EQ.0.and.ILEPTO.EQ.1)) Q2 = DBLE(Q2LP)
      PYVAR(33) =  0.0
      PYVAR(34) =  0.0

      IP=MAX(IPU1,IPU2)
      NS=N
c      write(6,*) ' in rgremn ',IPU1,IPU2
c      call pylist(1)
C...DEFINE INITIAL PARTONS, INCLUDING PRIMORDIAL KT
   10 DO 30  I=3,4
         IF(I.EQ.3) IPU=IPU1
         IF(I.EQ.4) IPU=IPU2
         K(I,1)=21
         K(I,3)=I-2
c         write(6,*) 'rgremn 1st :PARL(22),T2GKI ',PARL(22),T2gki,xfgki
         DO 20 J=1,5
   20    P(I,J)=0.D0
         IF(IPU.NE.0) THEN
            K(I,2)=K(IPU,2)
            P(I,5)=P(IPU,5)
c            P(I,5)=0.D0
c            write(6,*) 'rgremn: prkt1',prkt1,ptspl,phispl
            IF(IDIR.NE.0) THEN
               if(i.eq.3) then
                  if(IPY(40+I-2).eq.22) then
c for photon select gaussian,dkt^2/(kt0^2+kt^2) or dkt^2/(kt0^2+kt^2)^2
                     IF(IGAMKT.eq.1) Then
                        prkt = sngl(PRKT1)
                        CALL LPRIKT(prkt,PTSPL,PHISPL)
                     ELSEIF(IGAMKT.EQ.2) Then
                        pktmin=prkt1**2
                        pktmax=pktm**2+prkt1**2
                        PTSPL = sqrt(max(0d0, (pktmin*(pktmax/pktmin)**
     +                  draprn()-prkt1**2)))
                        PHISPL = 6.2832*draprn()
                     ELSEIF(IGAMKT.EQ.3) Then
                        pktmin=prkt1**2
                        pktmax=pktm**2+prkt1**2
                        PTSPL = sqrt(max(0d0, (pktmin*pktmax/(pktmax-
     +                  draprn()*pktm**2)-prkt1**2)))
                        PHISPL = 6.2832*draprn()
                     ELSE
                        write(6,*) ' rgremn no valid distr. for prim '
     +                  //'kt'
                        write(6,*) ' in photon selected: IGAMKT=',
     +                  IGAMKT
                        STOP
                     ENDIF

                  else
                     prkt = sngl(PRKT1)
                     CALL LPRIKT(prkt,PTSPL,PHISPL)
                  endif
               else
                  prkt = sngl(PRKT2)
                  CALL LPRIKT(prkt,PTSPL,PHISPL)
               endif
            ELSEIF(IDIR.EQ.0) THEN
               PTSPL=0.
               if(i.eq.3) then
                  if(IPY(40+I-2).eq.22) then
c for photon select gaussian,dkt^2/(kt0^2+kt^2) or dkt^2/(kt0^2+kt^2)^2
                     IF(IGAMKT.eq.1) Then
                        prkt = sngl(PRKT1)
                        CALL LPRIKT(prkt,PTSPL,PHISPL)
                     ELSEIF(IGAMKT.EQ.2) Then
                        pktmin=prkt1**2
                        pktmax=pktm**2+prkt1**2
                        PTSPL = sqrt(max(0d0, (pktmin*(pktmax/pktmin)**
     +                  draprn()-prkt1**2)))
                        PHISPL = 6.2832*draprn()
                     ELSEIF(IGAMKT.EQ.3) Then
                        pktmin=prkt1**2
                        pktmax=pktm**2+prkt1**2
                        PTSPL = sqrt(max(0d0, (pktmin*pktmax/(pktmax-
     +                  draprn()*pktm**2)-prkt1**2)))
                        PHISPL = 6.2832*draprn()
                     ELSE
                        write(6,*) ' rgremn no valid distr. for prim '
     +                  //'kt'
                        write(6,*) ' in photon selected: IGAMKT=',
     +                  IGAMKT
                        STOP
                     ENDIF
                  else
                     prkt = sngl(PRKT1)
c                     prkt = 0
                     CALL LPRIKT(prkt,PTSPL,PHISPL)
                  endif
               else
                  prkt = sngl(PRKT2)
                  CALL LPRIKT(prkt,PTSPL,PHISPL)
                  IF(IALMKT.EQ.1.AND.IPRO.EQ.12)  
     +            CALL RALMKT(SNGL(PARL(21) ),PTSPL,PHISPL)
               endif

            ENDIF
            IF(ISEMIH.EQ.1) THEN
               PTSPL=0.
               P(I,5)=0.D0
            ENDIF
            if(Ilha.ne.0)  then 
cccc                P(I,5)=0.D0
                PTSPL =0
            endif
            IF(I.EQ.3) GGPT = PTSPL
            IF(I.EQ.4) GPPT = PTSPL
c            write(6,*) 'rgremn: prkt1',prkt1,ptspl**2,phispl
c            write(6,*) 'rgremn: prkt2',prkt2
c            write(6,*) 'rgremn PT2gen,phigen ',pt2gen,phigen
            P(I,1)=DBLE(PTSPL*COS(PHISPL))
            P(I,2)=DBLE(PTSPL*SIN(PHISPL))
            PMS(I-2)=P(I,5)**2+P(I,1)**2+P(I,2)**2
c            write(6,*) 'rgremn: PMS=',I,PMS(I-2),P(I,5),P(I,1),P(I,2)
c        write(6,*) ' ILEP ',ILEP
         ELSE
            K(I,2)=K(IQ,2)
            IF(ILEPTO.EQ.0) THEN
               P(I,5)= 0.0D0
cccc               pms(i-2)=0.0D0

            ELSEIF(ILEPTO.EQ.1) THEN
               P(I,5)=-SQRT(Q2)
               PMS(I-2)= - Q2
c. shs = w**2 but calculated for photoprod.
               SHS= (- Q2 + PARL(22)+PYVAR(7-I)**2)
cccc         write(6,*) 'rgremn  :SHS,PYVAR(7-I)',SHS,PYVAR(7-I)
c              write(6,*) ' SHS ',SHS,' SHS_n =',SHSN
               IF(IDIR.EQ.0) THEN
                  SHS = DBLE(-Q2 + PARL(22) + DBLE(T2GKI))
c         write(6,*) 'rgremn new :SHS,PYVAR(7-I)',SHS,T2GKI
c         write(6,*) 'rgremn:PARL(22),T2GKI ',PARL(22),T2gki,xfgki
               ENDIF
            ENDIF
         ENDIF
   30 CONTINUE
c       write(6,*) ' pms(1),pms(2) ',pms(1),pms(2)
c      pause ' rgremn: after 30 '
c      write(6,*) 'before construction of partons ',ILEP
c      call pylist(1)
C...KINEMATICS CONSTRUCTION FOR INITIAL PARTONS
      IF(ILEP.EQ.0) then
         SHS=PYVAR(31)*PYVAR(32)*PYVAR(2)
     +       +(P(3,1)+P(4,1))**2+(P(3,2)+P(4,2))**2
      ENDIF
c      write(6,*) ' shs ,shat,t2gki ',shs,shat,t2gki
c      write(6,*) ' PYVAR 31,32,2 :',pyvar(31),PYVAR(32),PYVAR(2)
c      write(6,*) ' PMS(1),PMS(2) :',PMS(1),PMS(2)
c      write(6,*) ' P(3,1),P(3,2): ',P(3,1),P(3,2)
c      write(6,*) ' P(4,1),P(4,2): ',P(4,1),P(4,2)
c      write(6,*) ' Q2,PARL(22)',Q2,PARL(22)
c      SHS = SHAT
c      call pylist(1)
ccc      pause
      SHR=DSQRT(DMAX1(0.D0,SHS))
      IF(ILEP.EQ.0) THEN
         IF((SHS-PMS(1)-PMS(2))**2-4.D0*PMS(1)*PMS(2).LE.0.) GOTO 10
         P(3,4)=0.5D0*(SHR+(PMS(1)-PMS(2))/SHR)
         P(3,3)=DSQRT(DMAX1(0.D0,P(3,4)**2-PMS(1)))
         P(4,4)=SHR-P(3,4)
         P(4,3)=-P(3,3)
      ELSEIF(ILEP.EQ.1) THEN
         P(3,4)=P(IQ,4)
         P(3,3)=P(IQ,3)
         P(4,4)=P(IP,4)
         P(4,3)=P(IP,3)
c         write(6,*) ' rgremn IQ: ',P(IQ,3),P(IQ,4)
c         write(6,*) ' rgremn IP: ',P(IP,3),P(IP,4)
      ELSEIF(ILEP.EQ.2) THEN
         P(3,4)=P(IP,4)
         P(3,3)=P(IP,3)
         P(4,4)=P(IQ,4)
         P(4,3)=P(IQ,3)
      ENDIF
C...TRANSFORM PARTONS TO OVERALL CM-FRAME (NOT FOR LEPTOPRODUCTION)
      IF(ILEP.EQ.0) THEN
         MSTU(1)=3
         MSTU(2)=4
ccc         write(6,*) ' in ILEP=0 ',P(3,3),P(3,4),P(4,3),P(4,4)
c         call pylist(1)
         DROBO(3)=(P(3,1)+P(4,1))/SHR
         DROBO(4)=(P(3,2)+P(4,2))/SHR
         CALL pyrobo(MSTU(1),MSTU(2),0.D0,0.D0,-DROBO(3),-DROBO(4),0.D0)
c         CALL pylist(1)
         ROBO(2)=pyangl(P(3,1),P(3,2))
         CALL pyrobo(MSTU(1),MSTU(2),0.D0,-ROBO(2),0.D0,0.D0,0.D0)
c         CALL pylist(1)
         ROBO(1)=pyangl(P(3,3),P(3,1))
         CALL pyrobo(MSTU(1),MSTU(2),-ROBO(1),0.D0,0.D0,0.D0,0.D0)
c         CALL pylist(1)
ccc         write(6,*) 'IPY(47),IPU1,IPU2',IPY(47),IPU1,IPU2
         MSTU(2)=MAX(IPY(47),IPU1,IPU2)
         CALL pyrobo(MSTU(1),0,ROBO(1),ROBO(2),DROBO(3), DROBO(4),0D0)
         DROBO(5)=DBLE(MAX(-0.999999,MIN(0.999999,(PYVAR(31)-PYVAR(32))/
     +   (PYVAR(31)+PYVAR(32)))))
         CALL pyrobo(MSTU(1),0,0.D0,0.D0,0.D0,0.D0,DROBO(5))
c         CALL pylist(1)
c       ROBO(2) = pyangl(P(1,1),P(1,2))
c       ROBO(1) = pyangl(P(1,3),P(1,1))
         MSTU(1)=0
         MSTU(2)=0
      ENDIF
ccc      write(6,*) ' after boost : ILEP=',ILEP
ccc      CALL pylist(1)
ccc      pause 'after boost'
C...CHECK INVARIANT MASS OF REMNANT SYSTEM:
C...HADRONIC EVENTS OR LEPTOPRODUCTION
      IF(ILEP.LE.0) THEN
c here is changed; now according to PYTHIA 5.6
         PSYS(0,4)=P(3,4)+P(4,4)
         PSYS(0,3)=P(3,3)+P(4,3)
         PMS(0)=DMAX1(0.D0,PSYS(0,4)**2 - PSYS(0,3)**2)
         PMMIN0=DSQRT(PMS(0))
         PMMIN1=P(1,5)**2 + P(3,1)**2 + P(3,2)**2
         PMMIN1=DSQRT(PMMIN1)
         PMMIN2=P(2,5)**2 + P(4,1)**2 + P(4,2)**2
         PMMIN2=DSQRT(ABS(PMMIN2))
         PSYS(1,4) = DBLE(0.5*PYVAR(1)*(1. - PYVAR(30+1)))
         PSYS(2,4) = DBLE(0.5*PYVAR(1)*(1. - PYVAR(30+2)))
         IF(IRES(1).EQ.0.AND.(IABS(K(1,2)).EQ.11.OR.IABS(K(1,2)).EQ.13)
     +   ) PSYS(1,3)=PSYS(1,4)
         IF(PMMIN0+PMMIN1+PMMIN2.GT.PYVAR(1).OR.PMMIN1.GT.PSYS(1,4)
     +   .OR. PMMIN2.GT.PSYS(2,4)) THEN
            IPY(48)=1
ccc            Write(6,*) 'rgremn: IPY(48) ',PYVAR(1),PSYS(1,4),PSYS(2,4)
ccc            CALL PYLIST(1)
            RETURN
         ENDIF
      ELSE
         PEI=P(IQ,4)+P(IP,4)
         PZI=P(IQ,3)+P(IP,3)
         PMS(ILEP)=DMAX1(0.D0, PEI**2-PZI**2+P(5-ILEP,1)**2+
     +   P(5-ILEP,2)**2)
c         write(6,*) ' pms ',pms(ilep),q2
         IF(IPRO.EQ.12) THEN
c            write(6,*) ' PEI ',PEI,' PZI ',PZI
c            write(6,*) ' q2 ',Q2,dot1(IQ,IQ)
         ENDIF
         PMMIN=P(3-ILEP,5)+PYMASS(K(5-ILEP,2))+DSQRT(PMS(ILEP))
c         write(6,*) ' 1st PMMIN ',PMMIN
c         IF(IDIR.EQ.0) PMMIN=DBLE(T2GKI)
c     +    +PYMASS(K(5-ILEP,2))+DSQRT(PMS(ILEP))
c         write(6,*) ' 2st PMMIN ',PMMIN
c         write(6,*) ' PMS(ILEP) ',PMS(ILEP),PYMASS(K(5-ILEP,2)),T2GKI
c hju
         IF(IDIR.EQ.1) THEN
            PMMIN=PMMIN+PYPAR(12)
         ELSE
            PMMIN=DBLE(T2GKI)+PYMASS(K(5-ILEP,2))+DSQRT(PMS(ILEP)
     +      )
         ENDIF

chju         IF(SHR.LE.PMMIN+PYPAR(12)) THEN
         IF(SHR.LE.PMMIN) THEN
c             write(6,*) ' SHR < PMMIN '
c            WRITE(6,*) ' IN rgremn: IPY(48)=1 ',idir,ipro
c            WRITE(6,*) 'SHR  ',SHR
c            WRITE(6,*) 'PMS(ILEP)=',PMS(ILEP),' PMMIN=',PMMIN
c            WRITE(6,*) 'PYPAR(12)= ',PYPAR(12),'IQ=',IQ,' IP=',IP
c            WRITE(6,*) 'PYVAR(32 =',PYVAR(32)
ccc            pause 'rgremn: IPY(48) RETURN'
ccc            CALL pylist(1)
            IPY(48)=1
            RETURN
         ENDIF
      ENDIF
c      write(6,*) ' before subdivide remnant'
cc      CALL pylist(1)
C...SUBDIVIDE REMNANT IF NECESSARY, STORE FIRST PARTON
c      pause ' rgremn: bef 40 '
      NTRY = 0
   40 I=NS-1
      DO 70  JT=1,2
         ISN(JT)=0
c      write(6,*) ' KFL = ',KFL(1,JT)
         IF((ILEP.EQ.JT).OR. (IPY(40+JT)
     +   .EQ.22.AND.IRES(JT).EQ.0)) GOTO 70
cc         IF((IRES(JT).EQ.0.AND.ILEPTO.EQ.1).OR. (IPY(40+JT)
cc     +   .EQ.22.AND.IRES(JT).EQ.0)) GOTO 70
         IF(JT.EQ.1) IPU=IPU1
         IF(JT.EQ.2) IPU=IPU2
         CALL PYSPLI(IPY(40+JT),KFL(1,JT),KFLCH(JT),KFLSP(JT))
c       write(6,*) 'rgremn IPY=',IPY(40+JT),KFL(1,JT),KFLCH(JT),KFLSP(JT)
         I=I+2
         IS(JT)=I
         ISN(JT)=2
         K(I,1)=3
         K(I,2)=KFLSP(JT)
         K(I,3)=JT
         P(I,5)=PYMASS(K(I,2))
ccheat with masses for ILHA 
         if(ilha.gt.0) P(I,5)= DBLE(PYMASS( KFLCH(JT))) 
c         P(I,5) = 0
c         write(6,*) 'rgremn: I=',I,' mass ',P(I,5)
c          IF(IDIR.EQ.0) P(I,5) = 0.0
C...FIRST PARTON COLOUR CONNECTIONS AND TRANSVERSE MASS
         K(I+1,1)=-1
         K(I+1,3)=I
         K(I+1,2)=1000
         IF(IPY(34).GE.1) K(I+1,2)=1000+JT
         DO 50 J=1,5
   50    P(I+1,J)=0.D0
         IF(IPU.NE.0.AND.KFLSP(JT).EQ.21) THEN
            P(I+1,3)=IPU
            P(I+1,4)=IPU
            P(IPU+1,1)=I
            P(IPU+1,2)=I
            K(I,4)=IPU+IPU*MSTU(5)
            K(I,5)=IPU+IPU*MSTU(5)
            K(IPU,4)=MOD(K(IPU,4),MSTU(5))+I*MSTU(5)
            K(IPU,5)=MOD(K(IPU,5),MSTU(5))+I*MSTU(5)
         ELSEIF(IPU.NE.0) THEN
            IFLS=(3-ISIGN(1,KFLSP(JT)*(1102-IABS(KFLSP(JT)))))/2
            P(I+1,IFLS+2)=IPU
            P(IPU+1,3-IFLS)=I
            K(I,IFLS+3)=IPU
            K(IPU,6-IFLS)=MOD(K(IPU,6-IFLS),MSTU(5))+I*MSTU(5)
         ENDIF
         IF(KFLCH(JT).EQ.0) THEN
            P(I,1)=-P(JT+2,1)
            P(I,2)=-P(JT+2,2)
            IF(ISEMIH.EQ.1) THEN
               P(I,1) = - P(IP,1)
               P(I,2) = - P(IP,2)
            ENDIF
            PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2
cc          write(6,*) ' PMS = ',PMS(JT),JT,P(I,5)**2,P(I,1)**2+P(I,2)**2
         ELSE

C...WHEN EXTRA REMNANT PARTON OR HADRON: FIND RELATIVE PT, STORE
C...PRIMORDIAL KT SPLIT SHARED BETWEEN REMNANTS
c            write(6,*) ' look for prim. k_t '
            CALL LPRIKT(SNGL(PARL(14)),PTSPL,PHISPL)
c            write(6,*) ' PARL(14) =',PARL(14),PTSPL,PHISPL
            IF((IABS(K(JT,2)).EQ.11.OR.IABS(K(JT,2)).EQ.13).AND.
     +      IRES(JT).NE.0) THEN
cccccccccccccccc               CHI(JT)=(1. - yy)/(1. - PYVAR(30+JT))
               CHI(JT)=(1. - yy)/(1. - PYVAR(30+JT)*yy)
            ELSEIF((IABS(K(JT,2)).EQ.11.OR.IABS(K(JT,2)).EQ.13).AND.
     +      IRES(JT).EQ.0)  THEN
               CHI(JT)=YY
            ELSE
C...RELATIVE DISTRIBUTION OF ENERGY; EXTRA PARTON COLOUR CONNECTION
               If(IRspl.ge.1.and.IRspl.le.3) then
123               CALL lremh(0,KFLSP(JT),KFLCH(JT),CHI(JT))
                  icallz = icallz + 1
                  if(icallz.gt.30) then
                     write(6,*) ' rgremn icallz limit ',icallz
                     IPY(48)=1
                     RETURN
                  endif
                  if(chi(jt).gt.0.8.and.ilha.ne.0) goto 123
		   else
C...Using chosen fragmentation function in JETSET/PYTHIA
                  PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2
C	            write(6,*) ' rgremn pms ',pms(jt+4)  
			pmst=pms(jt+4)
                  icallz = 0
124               continue
                  icallz = icallz + 1
                  If(IRspl.eq.4) then 
			   CALL pyzdis(1,0,pmst,ZZ)
			elseif(IRspl.eq.5) then   
                     CALL pyzdis(1000,0,pmst,ZZ)
			else
			   write(6,*) ' rgremn IRspl not implemented ', IRspl
			   stop
			endif
                  if(icallz.gt.30) then
                     write(6,*) ' rgremn icallz limit '
                     IPY(48)=1
                     RETURN
                  endif
                  if(zz.gt.0.8.and.ilha.ne.0) goto 124
 			CHI(JT)=ZZ
		   endif
c               write(6,*) 'from PYremh:',KFLSP(JT),KFLCH(JT),CHI(JT)
            ENDIF
            P(I,1)=-P(JT+2,1)*(1.D0-DBLE(CHI(JT)))
     +       +DBLE(PTSPL*COS(PHISPL))
            P(I,2)=-P(JT+2,2)*(1.D0-DBLE(CHI(JT)))
     +       +DBLE(PTSPL*SIN(PHISPL))
            if(ilha.ne.0) then 
                 P(I,1)=0.
                 P(I,2)=0.
            endif
            IF(ISEMIH.EQ.1) THEN
               P(I,1)=-P(IP,1)*DBLE((1.-CHI(JT))+PTSPL*COS(PHISPL))
               P(I,2)=-P(IP,2)*DBLE((1.-CHI(JT))+PTSPL*SIN(PHISPL))
               if(ilha.ne.0) then 
                 P(I,1)=0.
                 P(I,2)=0.
               endif
            ENDIF

            PMS(JT+2)=P(I,5)**2+P(I,1)**2+P(I,2)**2
            I=I+2
            ISN(JT)=4
            DO 60 J=1,5
               K(I,J)=0
               K(I+1,J)=0
               P(I,J)=0.D0
   60       P(I+1,J)=0.D0
            K(I,1)=1
            K(I,2)=KFLCH(JT)
            K(I,3)=JT
            P(I,5)=PYMASS(K(I,2))
            P(I,1)=-P(JT+2,1)*DBLE(CHI(JT))-DBLE(PTSPL*COS(PHISPL))
            P(I,2)=-P(JT+2,2)*DBLE(CHI(JT))-DBLE(PTSPL*SIN(PHISPL))
            if(ilha.ne.0) then 
                  P(I,1)=0
                  P(I,2)=0
            endif
           IF(ISEMIH.EQ.1) THEN
               P(I,1)=-P(IP,1)*DBLE(CHI(JT)-PTSPL*COS(PHISPL))
               P(I,2)=-P(IP,2)*DBLE(CHI(JT)-PTSPL*SIN(PHISPL))
                if(ilha.ne.0) then 
                  P(I,1)=-P(IPU,1)
                  P(I,2)=-P(IPU,2)
               endif
            ENDIF
            PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2
C...end of update
            PMS(JT)=PMS(JT+4)/DBLE(CHI(JT))+PMS(JT+2)/
     +               (1.D0-DBLE(CHI(JT)))
c            write(6,*) ' jt,pms(jt),chi ',jt,pms(jt),chi(JT),yy
c            pause
            K(I+1,1)=-1
            K(I+1,3)=I
            K(I+1,2)=1000
            IF(IPY(34).GE.1) K(I+1,2)=1000+JT
            IF((IABS(KFLCH(JT)).GE.1.AND.IABS(KFLCH(JT)).LE.8).OR.
     +      IABS(KFLCH(JT)).EQ.21.OR.pycomp(IABS(KFLCH(JT))).EQ.90)
     +      THEN
               IFLS=(3-ISIGN(1,KFLCH(JT)*(1102-IABS(KFLCH(JT)))))/2
               P(I+1,IFLS+2)=IPU
               P(IPU+1,3-IFLS)=I
               K(I,1)=3
               K(I,IFLS+3)=IPU
               K(IPU,6-IFLS)=MOD(K(IPU,6-IFLS),MSTU(5))+I*MSTU(5)
            ELSE
               IF(IPY(34).GE.1) THEN
                  K(I,1)=1
                  K(I,3)=JT
               ENDIF
            ENDIF
         ENDIF
c           write(6,*) ' jt,pms(jt),pms(1) ',jt,pms(jt),pms(1)
   70 CONTINUE
c      write(6,*) ' ILEP ',ILEP
c      pause 'rgremn: before goto 40'
      IF(ILEP.GE.1.AND.
     +   SHR.LE.(DSQRT(PMS(1))+DSQRT(DABS(PMS(2))))) THEN
         NTRY = NTRY + 1
         IF(NTRY.LT.500) GOTO 40
c         write(6,*) 'rgremn IPY(48)=2 shr,dsqrt(pms(1)),dsqrt(pms(2))'
c     +    ,shr,dsqrt(pms(1)),dsqrt(pms(2))
c         write(6,*) ' rgremn pt2gen,shr/4 ',pt2gen,shr/4.
c         call lulist(1)
         IPY(48) = 2
         RETURN
      ENDIF


      IF(ILEP.EQ.0.AND.(PMS(1).GT.PSYS(1,4)**2.OR.
     +   PMS(2).GT.PSYS(2,4)**2)) THEN
         NTRY = NTRY + 1
         IF(NTRY.LT.500) GOTO 40
         IPY(48) = 2
         RETURN
      ENDIF
c      write(6,*) 'rgremn: after goto 40'
c      pause ' rgremn: bef kin '
c      call pylist(1)
c      pause
      N=I+1
C...RECONSTRUCT KINEMATICS OF REMNANTS
c      write(6,*) 'before remnant pms(1),pms(2) ',pms(1),pms(2)
      DO 80  JT=1,2
         IF((ILEPTO.EQ.1.AND.IRES(JT).EQ.0).OR.(IPY(40+JT)
     +   .EQ.22.AND.IRES(JT).EQ.0)) GOTO 80
         IF(ILEP.EQ.1) THEN
c            write(6,*) ' SHR = ',SHR,' PMS(JT)=',PMS(JT),PMS(3-JT)
c            write(6,*) ' Q2 = ',Q2,' T2 ',T2GKI
c            write(6,*) ' JT = ',JT
            PE=0.5D0*(SHR+(PMS(JT)-PMS(3-JT))/SHR)
            PZ=DSQRT(PE**2-PMS(JT))
c            write(6,*) ' PE,PZ ',PE,PZ
         ELSE
c new ala PYTHIA 5.6
            PE=PSYS(JT,4)
            PZ=DSQRT(PE**2-PMS(JT))
         ENDIF
         IF(KFLCH(JT).EQ.0) THEN
            P(IS(JT),4)=PE
            P(IS(JT),3)=PZ*(-1)**(JT-1)
         ELSE
            PW1=DBLE(CHI(JT))*(PE+PZ)
            P(IS(JT)+2,4)=0.5D0*(PW1+PMS(JT+4)/PW1)
            P(IS(JT)+2,3)=0.5D0*(PW1-PMS(JT+4)/PW1)*(-1)**(JT-1)
            P(IS(JT),4)=PE-P(IS(JT)+2,4)
            P(IS(JT),3)=PZ*(-1)**(JT-1)-P(IS(JT)+2,3)
         ENDIF
         PSYS(JT,3)=PZ*(-1)**(JT-1)
   80 CONTINUE


C...HADRONIC EVENTS: BOOST REMNANTS TO CORRECT LONGITUDINAL FRAME
      IF(ILEP.LE.0) THEN
c	write(6,*) ' before long boost '
c      CALL PYLIST(1)
ccc      pause
C new a la PYTHIA 5.4.....
C... check if longitudinal boost needed - if so pick two systems.
         PDEV = ABS(PSYS(0,4)+PSYS(1,4)+PSYS(2,4) -
     +   PYVAR(1) + ABS(PSYS(0,3)+PSYS(1,3)+PSYS(2,3)))

         IF(PDEV.LE.1.E-6*PYVAR(1).OR. (K(1,2).EQ.22.AND.IRES(1).EQ.0))
     +   RETURN
c           write(6,*) ' now do the boost ',PDEV
c           CALL pylist(1)
         IF(ISN(1).EQ.0) THEN
            IR=0
            IL=2
         ELSEIF(ISN(2).EQ.0) THEN
            IR=1
            IL=0
         ELSEIF((1.-PYVAR(31)).GT.0.2.AND.(1.-PYVAR(32)).GT.0.2) THEN
            IR=1
            IL=2
         ELSEIF((1.-PYVAR(31)).GT.0.2) THEN
            IR=1
            IL=0
         ELSEIF((1.-PYVAR(32)).GT.0.2) THEN
            IR=0
            IL=2
         ELSEIF(PMS(1)/PSYS(1,4)**2.GT.PMS(2)/PSYS(2,4)**2) THEN
            IR=1
            IL=0
         ELSE
            IR=0
            IL=2
         ENDIF
         IG=3-IR-IL

C...Construct longitudinal boost
         IF((IG.EQ.1.AND.ISN(1).EQ.0).OR.(IG.EQ.2.AND.ISN(2).EQ.0))
     +   THEN
            PPB=PYVAR(1)
            PNB=PYVAR(1)
         ELSE
            PPB=PYVAR(1)-PSYS(IG,4)-PSYS(IG,3)
            PNB=PYVAR(1)-PSYS(IG,4)+PSYS(IG,3)
         ENDIF

         PMTB=PPB*PNB
         PMTR=PMS(IR)
         PMTL=PMS(IL)
         SQLAM=DSQRT(DMAX1(0.D0,(PMTB-PMTR-PMTL)**2-4.D0*PMTR*PMTL))
         SQSGN=DSIGN(1.D0,PSYS(IR,3)*PSYS(IL,4)-PSYS(IL,3)*PSYS(IR,4))
         RKR=(PMTB+PMTR-PMTL+SQLAM*SQSGN)/(2.D0*(PSYS(IR,4)+PSYS(IR,3))*
     +   PNB)
         RKL=(PMTB+PMTL-PMTR+SQLAM*SQSGN)/(2.D0*(PSYS(IL,4)-PSYS(IL,3))*
     +   PPB)
         BER=(RKR**2-1.D0)/(RKR**2+1.D0)
         BEL=-(RKL**2-1.D0)/(RKL**2+1.D0)

C...Perform longitudinal boost
         MSTU(1)=0
         MSTU(2)=0
         IF(IR.EQ.1) THEN
            CALL pyrobo(IS(1),IS(1)+ISN(1)-2,0.D0,0.D0,0D0,0D0,BER)
         ELSE
            CALL pyrobo(3,NS,0.D0,0.D0,0D0,0D0,DBLE(BER))
         ENDIF
         IF(IL.EQ.2) THEN
            CALL pyrobo(IS(2),IS(2)+ISN(2)-2,0.D0,0.D0,0D0,0D0,BEL)
         ELSE
            CALL pyrobo(3,NS,0.D0,0.D0,0D0,0D0,BEL)
         ENDIF
ccc         call pylist(1)
ccc         pause ' rgremn pythia end '
C...LEPTOPRODUCTION EVENTS: BOOST COLLIDING SUBSYSTEM
      ELSE

         IMIN=21
         MSTU(1) = IMIN
         IMAX=MAX(IP,IPY(47))
         MSTU(2) = IMAX
c         write(6,*) ' rgremn: before long frame '
c         CALL pylist(1)
         PZF=PZ*(-1)**(ILEP-1)
         PT2=P(5-ILEP,1)**2+P(5-ILEP,2)**2
         PHIPT=pyangl(P(5-ILEP,1),P(5-ILEP,2))
         CALL pyrobo(IMIN,IMAX,0.D0,-PHIPT,0.D0,0.D0,0.D0)
c         write(6,*) ' rgremn: after phi rot '
c         CALL pylist(1)
         RQP=P(IQ,3)*(PT2+PEI**2)-P(IQ,4)*PEI*PZI
         SINTH=P(IQ,4)*DSQRT(PT2*(PT2+PEI**2)/(RQP**2+PT2* P(IQ,4)**2*
     +   PZI**2))*DSIGN(1.D0,-RQP)
         CALL pyrobo(IMIN,IMAX,DASIN(SINTH),0.D0,0.D0,0.D0,0.D0)
c         write(6,*) ' rgremn: after theta rot '
c         CALL pylist(1)
         DBETAX=(-PEI*PZI*SINTH+ DSQRT(PT2*(PT2+PEI**2-(PZI*SINTH)*
     +   *2)))/ (PT2+PEI**2)
         CALL pyrobo(IMIN,IMAX,0.D0,0.D0,DBETAX,0.D0,0.D0)
         CALL pyrobo(IMIN,IMAX,0.D0,PHIPT,0.D0,0.D0,0.D0)
c         write(6,*) ' rgremn: after x boost and phi '
c         CALL pylist(1)
         PEM=P(IQ,4)+P(IP,4)
         PZM=P(IQ,3)+P(IP,3)
c            write(6,*) ' PEM ',PEM,' PZM ',PZM
         DBETAZ=(-PEM*PZM+
     +   PZF*DSQRT(PZF**2+PEM**2-PZM**2))/(PZF**2+PEM**2)
c         write(6,*) ' boost check DBETAZ=',DBETAZ
c         write(6,*) 'PZM,PZF,PEM ',PZM,PZF,PEM
         CALL pyrobo(IMIN,IMAX,0.D0,0.D0,0.D0,0.D0,DBETAZ)
c         write(6,*) ' rgremn: after DBETZ '
c         CALL pylist(1)
         MSTU(1) = 0
         MSTU(2) = 0
         CALL pyrobo(3,4,DASIN(SINTH),0.D0,DBETAX,0.D0,0.D0)
         CALL pyrobo(3,4,0.D0,PHIPT,0.D0,0.D0,DBETAZ)
         PEM=P(IQ,4)+P(IP,4)
         PZM=P(IQ,3)+P(IP,3)
c            write(6,*) ' now after last boost '
c            write(6,*) ' PEM ',PEM,' PZM ',PZM
c         write(6,*) ' rgremn: end '
c         CALL pylist(1)
c         pause 'rgremn: END '
c         CALL pylist(1)
      ENDIF
c         write(6,*) ' rgremn: end '
c      call pylist(1)
      RETURN
      END
