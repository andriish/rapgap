      SUBROUTINE ppmeps
      IMPLICIT NONE
#include "rglujets.inc"
#include "rgrapgki.inc"
#include "rgpara1.inc"
#include "rgparas.inc"
#include "rgpart.inc"
#include "rgluco.inc"
#include "rgdisdif.inc"
#include "rgdiffr.inc"
#include "rgludat1.inc"
#include "rgludat3.inc"
#include "rgherup.inc"
#include "rgpypars.inc"
#include "rgcolcon.inc"
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
      DOUBLE PRECISION dboost
      REAL CUT,X,Y,W2,Q2LP,U
      DOUBLE PRECISION PARL
      Integer LLST
      COMMON /RAPTOU/ CUT(14),LLST(40),PARL(30),X,Y,W2,Q2LP,U
      DOUBLE PRECISION DBETA(2,3),DBETAR(2,3)
      DOUBLE PRECISION STHETA(2),SPHI(2),STHETAR(2),SPHIR(2)
      Double Precision PYPAR,PYVAR
      Integer IPY
      COMMON /RGPYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      REAL XPY,SH,TH,UH,Q2PY
      Integer ISUB,KFL
      COMMON /PYPROC/ ISUB,KFL(3,2),XPY(2),SH,TH,UH,Q2PY
      REAL XQPY
      COMMON /MYINT1/ XQPY(2,-6:6)
      REAL XPQ1,XPQ2
      DIMENSION XPQ1(-6:6),XPQ2(-6:6)
      DOUBLE PRECISION PS(40,5),PSAVE(2,5)
	Integer KPOM(2)
      REAL XPRO
	DOUBLE PRECISION QMAX
      Integer LST,IRES
      COMMON/EPPARA/LST(30),IRES(2)
      Integer KS,NS
      COMMON /COLR/ KS(40,5),NS
      REAL SNGL

      Double Precision remPARJ32
      Double Precision remPARJ(5)
      Integer I,J,NCALL,NPRIN,NCHECK,NIPH,NB2,NSTB,IFL,IT,IPU1,IPU2
      Integer IRFL,IPOM,MSTJO,IPROT
      Double Precision SHRAP,Q2LM,UMAS2,XR,DETOTR,XGTEST,Q12
      Double Precision PC1,PC2,PC3,PC4,CHEC,P1,P2,P3,P4
      Double Precision RMX2
	Integer IBEAM(2)
      
      Double Precision cmspart(5),sh_test
      Double Precision DOT
	Integer MSTJO14
      Integer kpom1,kpom2,npo
      Double Precision MT1,MT2
      integer IFINAL
      Integer MINT
      Double Precision VINT
      COMMON/PYINT1/MINT(400),VINT(400)
      
      logical first,LPOM
      data    first /.TRUE./
c

      DATA IFINAL/1/
      DATA NCALL/0/
      DATA NPRIN/0/
      DATA NCHECK/0/
      If(first) Then
         first = .FALSE.
         remPARJ32 = PARJ(32)
      ELSE
         PARJ(32)=remPARJ32
      EndIf
      DO I=1,40
         DO J=1,5
            KS(I,J) = 0
            PS(I,J) = 0.D0
         ENDDO
      ENDDO

      NCALL = NCALL + 1
      IPY(48)=0
      LST(21)=0
	
C ILEPTO=1 parton shower a la LEPTO
C ILEPTO=0 parton shower a la PYTHIA
      ILEPTO = 0
      SHRAP=DBLE(SHH)
      Q2PY = SNGL(Q2Q)
      Q2LP = SNGL(Q2)
      Q2LM = Q2
c change to e in + z direction
      DO 20  I=1,N
         P(I,3) = - P(I,3)
   20 continue
      Do I=N+1,100
	Do J=1,5
         P(I,J) = 0
	   K(I,J) = 0
	ENDDO
	ENDDO
	
      UMAS2 = PYMASS(K(2,2))**2
      NIPH = NIA1
      NB2 = 2
      IF(IRES(1).EQ.1) THEN
         NB2  = NIA2
      ELSE
         NIPH = NIA1
      ENDIF
      IF(IDIR.EQ.0) THEN
         NB2=NIA1+1
         UMAS2 = DBLE(T2GKI)
      ENDIF
      PARL(21) = 2.D0*DOT1(1,2)
      IF(IDIR.EQ.0) then
          PARL(21) = PARL(21) * DBLE(XFGKI)
      ENDIF
      PARL(21) = PARL(21) + PYMASS(K(1,2))**2 +UMAS2
c.hju parl(22) = parton - parton CMS
      PARL(22) = 2.D0*DOT1(NB2,NIPH)

      IF(IABS(K(1,2)).EQ.22.OR.IABS(K(1,2)).EQ.23) Q2LM=0.0D0

      IF(NCALL.LE.NPRIN) CALL pylist(1)
      IF(NCHECK.EQ.1) THEN
         CALL pylist(1)
      ENDIF



	IBEAM(1)=1
	IBEAM(2)=2
      NPO = 0
      IF(IDIR.EQ.0) THEN
c find pomeron
         DO I=1,N
	   IF(K(I,2).EQ.990.OR.IABS(K(I,2)).EQ.211.OR.
     +	   IABS(K(I,2)).EQ.210) THEN
            IF(K(I,1).EQ.21) THEN
		IF(K(I,3).EQ.1) IBEAM(1)=I
		IF(K(I,3).EQ.2) IBEAM(2)=I
		ENDIF
         ENDIF
	   ENDDO
	   nstb = 0
         DETOTR=P(IBEAM(1),4)+P(IBEAM(2),4)
         DBETAR(2,1)=(P(IBEAM(1),1)+P(IBEAM(2),1))/DETOTR
         DBETAR(2,2)=(P(IBEAM(1),2)+P(IBEAM(2),2))/DETOTR
         DBETAR(2,3)=(P(IBEAM(1),3)+P(IBEAM(2),3))/DETOTR
         CALL pyrobo(nstb,0,0.D0,0.D0,-DBETAR(2,1),-DBETAR(2,2),
     +   -DBETAR(2,3))
         SPHIR(2)=pyangl(P(IBEAM(1),1),P(IBEAM(1),2))
         CALL pyrobo(nstb,0,0.D0,-SPHIR(2),0.D0,0.D0,0.D0)
         STHETAR(2)=pyangl(P(IBEAM(1),3),P(IBEAM(1),1))
         CALL pyrobo(nstb,0,-STHETAR(2),0.D0,0.D0,0.D0,0.D0)
         XPRO = XPR/XFGKI
	   NPO=0
	   IF(IBEAM(1).EQ.1) THEN
		KPOM(1) = 0
		K(3,1)=0
	   ELSEIF(IBEAM(1).NE.1) THEN
		NPO = NPO+1
		KPOM(1)=1
	      DO I=1,5
		KS(15,I) = K(IBEAM(1),I)
	      PS(15,I) = P(IBEAM(1),I)
		KS(16,I) = K(3,I)
	      PS(16,I) = P(3,I)
	      END DO
		K(3,1)=0
		ENDIF
	      IF(IBEAM(2).EQ.2) THEN
		KPOM(2) = 0
	      ELSEIF(IBEAM(2).NE.2) THEN
		NPO = NPO+1
		KPOM(2) = 1
	      DO I=1,5
		KS(17,I) = K(IBEAM(2),I)
	      PS(17,I) = P(IBEAM(2),I)
		KS(18,I) = K(N,I)
	      PS(18,I) = P(N,I)
	      END DO
		K(N,1)=0
		ENDIF
      ELSE
         XPRO = XPR
      ENDIF
c correction necessary for larger t, recalcualte xpro assuming collinear partons
c included 26.12. 2009
      IF(IDIR.EQ.0) then
         do i=1,4
           cmspart(i)=p(nia1,i)+p(nia2,i)
         end do
         sh_test= dot(cmspart,cmspart)
         xpro = sh_test/4./xel/P(1,4)/P(2,4)/xfgki
      endif
      PYVAR(32) = XPRO
      XR = DBLE(XPRO)
      PYVAR(31) = XEL
      DO I=1,5
         PSAVE(1,I)=P(1,I)
         PSAVE(2,I)=P(2,I)
      ENDDO
      DBETA(2,1)=0d0
      DBETA(2,2)=0d0
      SPHI(2)=0d0
      STHETA(2)=0d0
      dboost = (PYVAR(31)-PYVAR(32)) /(PYVAR(31)+PYVAR(32))
      DBETA(2,3)=dboost
C...Transform to gamma g cms, boost parameters in double precision.
      CALL pyrobo(3,0,0d0,0d0,0.0d0,0.0d0,-dboost)
      IF(IDIR.EQ.0) THEN
         DO J=1,5
            KS(31,J) = K(NIA1+1,J)
            PS(31,J) = P(NIA1+1,J)
            KS(32,J) = K(N,J)
            PS(32,J) = P(N,J)
         ENDDO
      ENDIF

      IF(NCALL.LE.NPRIN) CALL pylist(1)
      
C...Save event record in  cms
      DO 30 I=1,N
         DO 30 J=1,5
            KS(I,J)=K(I,J)
   30 PS(I,J)=P(I,J)
C...Rearrange event record to rasspa standard
      DO 40 J=1,5
         K(5,J)=KS(NIA1,J)
         P(5,J)=PS(NIA1,J)
         K(6,J)=KS(NIA2,J)
         P(6,J)=PS(NIA2,J)
         K(8,J)=KS(NF1,J)
         P(8,J)=PS(NF1,J)
c
         K(11,J) = 0
         K(12,J) = 0

         K(10,J)=KS(NF2,J)
   40 P(10,J)=PS(NF2,J)
      K(5,3)=3
      K(6,3)=4
      K(7,3)=5
      K(8,3)=6
      K(9,3)=5
      K(10,3)=6
      DO 50 I=5,10
   50 K(I,1)=21
      K(9,1)=0
c end of redefintion
      K(6,2)=KS(NIA2,2)
      N=19
      
      
      IF(NCALL.LE.NPRIN) THEN
         CALL pylist(1)
      ENDIF
C...Partons with colour information in hadronic cms frame.
      DO 80  I=19,MSTU(4)
         DO 80 J=1,5
            K(I,J)=0
            P(I,J)=0.D0
   80 V(I,J)=0.
      NS=20

      DO 90  J=1,5
         N=NS+7
         K(NS+1,J)=K(5,J)
         P(NS+1,J)=P(5,J)
         K(NS+3,J)=K(6,J)
         P(NS+3,J)=P(6,J)
         K(NS+5,J)=K(8,J)
         P(NS+5,J)=P(8,J)
         K(NS+6,J)=K(10,J)
         P(NS+6,J)=P(10,J)
         IF(IRES(1).NE.0) THEN
            K(NS+1,1) = 13
         ENDIF
   90 CONTINUE
C...Old standard continuation lines
      K(NS+2,1)=-1
      K(NS+2,3)=NS+1
      K(NS+4,1)=-1
      K(NS+4,3)=NS+3
      P(NS+4,3)=27
      P(NS+4,4)=27
C...Origin and Colour info for incoming parton
      K(NS+3,1)=13
      K(NS+3,3)=2
      K(NS+3,4)=27
      K(NS+3,5)=27
C...Colour info for two outgoing partons
      K(NS+5,1)=3
      K(NS+6,1)=3
	If(ipro.eq.102.or.ipro.eq.500) K(NS+5,1)=1
C...Effective outgoing parton = sum of both outgoing partons
      K(NS+7,1)=14
      K(NS+7,3)=3

      CALL COLRES

      DO 100 J=1,4
         P(NS+7,J)=P(8,J)+P(10,J)
  100 CONTINUE
      P(NS+7,5)=DSQRT(DMAX1(0.D0, P(NS+7,4)**2-P(NS+7,1)**2-P(NS+
     +7,2)**2- P(NS+7,3)**2))
      N=NS+7

C...Scale for bremsstrahlung etc.

      SH = SNGL(SHRAP)
      IPY(40)=10
      IPY(47)=N
C...Save quantities for later use.
      IPY(41)=KINT(1,2)
      IPY(42)=KINT(2,2)
      PYVAR(2)=SNGL(PARL(21))
      PYVAR(3)=SNGL(P(1,5))
      PYVAR(4)=SNGL(P(2,5))
      PYVAR(1)=SQRT(PYVAR(2))
      PYVAR(5)=PYVAR(1)/2.
      XPY(1)=PYVAR(31)
      XPY(2)=PYVAR(32)
         IF(XPY(1).GT.1.1) THEN
            write(6,*) ' ppmeps : XPY(1) > 1. ',XPY(1)
            write(6,*) ' ppmeps :s_hat, W**2, x_p',
     +        SHRAP,PYVAR(2),PYVAR(32)
            write(6,*) ' ppmeps :y, sss, q2lp, xel, xel/yy ',
     +        yy,sss,q2lp,xel,xel/yy
         ENDIF
         XPY(1)=AMIN1(1.0,XPY(1))
         CALL RGSTFU(IPY(41),XPY(1),Q2PY,XPQ1)
         IF(NCALL.LE.NPRIN) THEN
            write(6,*) ' IPY(41),XPY(1),Q2PY,XPQ1', IPY(41),XPY(1),
     +      Q2PY,XPQ1
         ENDIF

c.hju CALL RGSTFU(K(2,2),XR,Q2,XPQ)
      IF(XPY(2).GT.1.) write(6,*) 'ppmeps : XPY(2) > 1. ',XPY(2)
      CALL RGSTFU(IPY(42),XPY(2),Q2PY,XPQ2)
      DO 110 I=-6,6
         IF(XPQ2(I).EQ.0.0) XPQ2(I)=1.E-11
  110 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' IPY(42),XPY(2),Q2PY,XPQ2',
     +    IPY(42),XPY(2),Q2PY,XPQ2
      ENDIF
      DO 120 IFL=-6,6
         XQPY(1,IFL)=XPQ1(IFL)
  120 XQPY(2,IFL)=XPQ2(IFL)
      ISUB=39
      IPY(11)=1
      KFL(2,1)=KS(NIA1,2)
      KFL(2,2)=K(6,2)
      KFL(1,1)=KFL(2,1)
      KFL(1,2)=KFL(2,2)
      KFL(3,1)=K(1,2)
      KFL(3,2)=K(27,2)


C...Generate timelike parton shower (if required)
      IF(IPY(13).EQ.1) THEN
         IF(IFINAL.EQ.1) THEN
            MT1 = P(25,1)**2 + P(25,2)**2 + P(25,5)**2
            MT2 = P(26,1)**2 + P(26,2)**2 + P(26,5)**2
            QMAX = (MT1 + MT2)/2.d0
            QMAX = 4d0*QMAX
         ELSEIF(IFINAL.EQ.2) THEN
            QMAX=P(27,5)**2
         ELSEIF(IFINAL.EQ.3) THEN
            MT1 = P(25,5)**2
            MT2 = P(26,5)**2
            QMAX = (MT1 + MT2)/2.d0
            QMAX = 4d0*QMAX
         ENDIF 
         QMAX = QMAX *SCALFA
         QMAX=DSQRT(QMAX)
         CALL pyshow(25,26,QMAX)
      ENDIF
      MINT(51) = 0
      MSTP(71) = IPY(13)
c check if there are resonances to decay
      if((K(25,2).gt.22.and.K(25,2).le.25)) Then
        CALL PYRESD(25)
      Endif
      if((K(26,2).gt.22.and.K(26,2).le.25)) then 
        CALL PYRESD(26)  
      Endif
      IF(MINT(51).EQ.1) then
             write(6,*) ' caps: problem with res decay '
      endif
      IT=25
      IF(N.GE.27) IT=27
      NS=N
C...Generate spacelike parton shower (if required)
      IPU1=0
changed
      IF(IRES(1).EQ.1) IPU1=21
      IPU2=23

      IF(IPY(14).GE.1) THEN
         CALL rasspa(IPU1,IPU2)
         IF(LST(21).NE.0) THEN
            IF(LST(21).EQ.55) THEN
            ELSE
               RETURN
            ENDIF
         ENDIF
      ENDIF
      IF(.NOT.(IPY(14).GE.1).OR.LST(21).EQ.55) THEN
         IF(LST(21).EQ.55) LST(21)=0
         DO 130 I=NS+1,NS+4
            DO 130 J=1,5
               K(I,J)=0
               P(I,J)=0.D0
  130    V(I,J)=0.
Cadded
         IF(IRES(1).EQ.0) THEN
            K(NS+1,1)=11
            K(NS+1,3)=21
         ELSE
            K(NS+1,1)=13
            K(NS+1,3)=21
            K(NS+1,4)=21
            K(NS+1,5)=21
            K(21,4)=K(21,4) + (NS+1)*MSTU(5)
            K(21,5)=K(21,5) + (NS+1)*MSTU(5)
         ENDIF
         K(NS+1,2)=KFL(2,1)
         DO 140 J=1,5
  140    P(NS+1,J)=P(21,J)
         K(NS+2,1)=-1
         K(NS+2,3)=NS+1
         K(NS+3,1)=13
         K(NS+3,2)=KFL(2,2)
         K(NS+3,3)=23
         K(NS+3,4)=23
         K(NS+3,5)=23
         P(NS+3,1)= 0.D0
         P(NS+3,2)= 0.D0
         P(NS+3,5)= 0.D0
         IF(ISEMIH.EQ.1) THEN
            Q12 = PS(NIA2,5)**2
            P(NS+3,1) = PS(NIA2,1)
            P(NS+3,2) = PS(NIA2,2)
            P(NS+3,5) = PS(NIA2,5)
         ELSE
            Q12 = 0.D0
         ENDIF
         IF(IRES(1).EQ.0) THEN
            IF(Q2LM.GT.2.D0) THEN
               P(NS+3,3)=(P(IT,5)**2+Q2LM+Q12)*(P(21,4)-P(21,3))/(2.D0*
     +         Q2LM)
               P(NS+3,4)=-P(NS+3,3)
            ELSE
               P(NS+3,3)=-(P(IT,5)**2+Q2LM+Q12)/(P(21,4)+P(21,3))/2.D0
               P(NS+3,4)=-P(NS+3,3)
            ENDIF
         ELSE
            P(NS+3,3) = P(23,3)
            P(NS+3,4) = P(23,4)
         ENDIF
         IF(NCALL.LT.NPRIN) THEN
            write(6,*) ' ppmeps NS+3,23,it ',P(NS+3,3),P(23,3),P(IT,5)
         ENDIF

         K(NS+4,1)=-1
         K(NS+4,3)=NS+3
         P(NS+4,3)=23
         P(NS+4,4)=23
         P(24,1)=NS+3
         P(24,2)=NS+3
         K(23,4)=K(23,4)+(NS+3)*MSTU(5)
         K(23,5)=K(23,5)+(NS+3)*MSTU(5)
         IPU1=0
         IF(IRES(1).EQ.1) IPU1=NS+1
         IPU2=NS+3
         N=N+4
      ENDIF
C...Rotate and boost outgoing parton shower

      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' after rotate and boost outgoing PS '
         call pylist(1)
      ENDIF
C...Hadron remnant and primordial kt
      IPY(47)=N
      IF(N.GT.MSTU(4)-20) THEN
         WRITE(6,*) ' ppmeps before rgremn: no more memory in PYJETS'
         LST(21)=51
         RETURN
      ENDIF
      CALL rgremn(IPU1,IPU2)
      IF(NCALL.LE.NPRIN) CALL pylist(2)
      IF(IPY(48).GE.1) THEN
         WRITE(6,*) ' ppmeps after rgremn IPY(48) = ',IPY(48)
         LST(21)=47+IPY(48)
         RETURN
      ENDIF
C...Rearrange partons along strings
      MSTU(24)=0
      MSTU(28)=0
      IF(IFPS.EQ.10) THEN
         MSTJO = MSTJ(105)
         MSTJ(105) = -1
      ENDIF
      IF(N.GT.MSTU(4)-20) THEN
         WRITE(6,*) ' ppmeps before pyprep: no more memory in PYJETS'
         LST(21)=51
         RETURN
      ENDIF
      IF(ILHA.GE.1) then
	   MSTJO14=MSTJ(14)
	   MSTJ(14)=-1   ! to avoid pyprep producing clusters
	ENDIF
cGB..
*==> Don't allow certain decay modes if collapsing
      RMX2     = -Q2 + PARL(22) + DBLE(T2GKI)
      LPOM=.NOT.(NG.EQ.20.OR.NG.EQ.21)
      
*
      If (RMX2.LT.4.0.AND.IDIR.EQ.0.AND.LPOM) Then
         If(ABS(KPA).LE.2) Then
            PARJ(32) = 0.35
            GOTO 180
         ElseIf(ABS(KPA).EQ.3) Then
            PARJ(32) = 0.1
            GOTO 180
         EndIf
  180    Continue
*
         remPARJ(1) = PARJ(11)
         remPARJ(2) = PARJ(12)
         remPARJ(3) = PARJ(13)
         remPARJ(4) = PARJ(27)
         remPARJ(5) = PARJ(28)
*
         PARJ(11) = 1.0
         PARJ(12) = 1.0
         PARJ(13) = 1.0
         PARJ(27) = 0.1
         PARJ(28) = 0.9
*
         CALL pyprep(0)
         PARJ(11) = remPARJ(1)
         PARJ(12) = remPARJ(2)
         PARJ(13) = remPARJ(3)
         PARJ(27) = remPARJ(4)
         PARJ(28) = remPARJ(5)
      elseif (RMX2.GT.9.0.AND.RMX2.LT.14.AND.
     + ABS(KPA).EQ.4.AND.IDIR.EQ.0) Then
         remPARJ(2) = PARJ(12)
         remPARJ(3) = PARJ(13)
         PARJ(13) = 1.0
         CALL pyprep(0)
         PARJ(13) = remPARJ(3)
         PARJ(12) = remPARJ(2)

      Else
         CALL pyprep(0)

      EndIf
      IF(ILHA.GE.1) MSTJ(14)= MSTJO14  
cGB....
      IF(IFPS.EQ.10) THEN
         MSTJ(105) = MSTJO
      ENDIF
      IF((MSTU(24).NE.0.OR.MSTU(28).NE.0.AND.ILHA.EQ.0).OR.
     &(MSTU(24).NE.0)) THEN
         WRITE(6,*) ' ppmeps: PYPREP error MSTU(24)= ',MSTU(24),MSTU(28)
         write(6,*) ' color configuration ICOLORA,IRESPRO ',
     +     ICOLORA,IRESPRO
         write(6,*) ' IRES(1)  ',IRES(1)
         call pylist(2)
         LST(21)=50
         RETURN
      ENDIF

C...Clean up event record -> order:
      IF(NCALL.LT.NPRIN) then
	   write(6,*) ' before clean up event '
	   CALL pylist(1)
	   endif
      DO 190 J=1,5
         K(N+2,J)=K(3,J)
         P(N+2,J)=P(3,J)
         K(N+1,J)=K(4,J)
         P(N+1,J)=P(4,J)
         K(10,J)=0
  190 CONTINUE
      DO 200 I=1,12
         DO 200 J=1,5
            PS(I,J) = P(I,J)
  200 KS(I,J) = K(I,J)
	IPOM = 0
	IF(NPO.GT.0.AND.KPOM(1).NE.0) THEN
        Ipom = ipom + 1
        kpom1 = 2 + ipom
      end if
	IF(NPO.GT.0.AND.KPOM(2).NE.0) THEN
        Ipom = ipom + 1
        kpom2 = 2 + ipom
      end if
      
      Nia1=3+ipom
      Nia2=5+ipom
      NF1=7+ipom
      NF2=8+Ipom
	K(25,3) = Nia1+1
	K(26,3) = Nia2+1
	do J=1,5
         K(nia1,J)=K(N+2,J)
         P(nia1,J)=P(N+2,J)
         K(Nia1+1,J)=K(21,J)
         P(Nia1+1,J)=P(21,J)
         K(nia2,J)=K(N+1,J)
         P(nia2,J)=P(N+1,J)
         K(Nia2+1,J)=K(23,J)
         P(Nia2+1,J)=P(23,J)
         K(NF1,J)=K(25,J)
         P(NF1,J)=P(25,J)
         K(NF2,J)=K(26,J)
         P(NF2,J)=P(26,J)	   
         IF(NPO.GT.0.AND.KPOM(1).NE.0) THEN
	   K(kpom1,J)=KS(15,J)
	   P(kpom1,J)=PS(15,J)
	   K(nf2+kpom1-2,J)=KS(16,J)
	   P(nf2+kpom1-2,J)=PS(16,J)
	   endif
	   IF(NPO.GT.0.AND.KPOM(2).NE.0) THEN
	   K(kpom2,J)=KS(17,J)
	   P(kpom2,J)=PS(17,J)
	   K(nf2+kpom2-2,J)=KS(18,J)
	   P(nf2+kpom2-2,J)=PS(18,J)
	   endif
		
	end do
	   K(3,1)=21
	   K(3,3)=1
	   K(4,1)=21
	   K(4,3)=3
	   K(5,1)=21
	   K(5,3)=2
	   K(6,1)=21
	   K(6,3)=5
         k(nf1,1)=21  !!!!!
         k(nf2,1)=21  !!!!!
	   IF(ipro.eq.102.or.ipro.eq.500) then 
            k(nf2,1)=-1
            nf2=-999
	      K(26,1)=0
	   Endif
      DO J=1,5
         P(1,J) = PSAVE(1,J)
         P(2,J) = PSAVE(2,J)
      end do
C...Zero irrelevant lines with K(I,1)<0
      DO 230 I=1,N
         IF(K(I,1).LT.0) THEN
            DO 220 J=1,5
               K(I,J)=0
  220       P(I,J)=0.D0
         ENDIF
  230 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' ppmeps: before PYEDIT '
         CALL pylist(1)
      ENDIF
      
C...Delete internal parton lines, i.e. with K(I,1)=13,14
      CALL pyedit(14)
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' ppmeps:after PYEDIT(14) '
         CALL pylist(1)
      ENDIF
C...Delete empty lines
      CALL pyedit(12)
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' ppmeps:after PYEDIT(12) '
         CALL pylist(1)
      ENDIF
      

c check for enegy momentum conservation
      DO 250 I=1,N
         IF(P(I,5).LT.0.0D0.OR.K(I,1).GT.2) GOTO 250
         DO 240 J=1,4
            IF(ABS(P(I,J)).LE.1.E-6) THEN
               P(I,J)=0.0D0
            ENDIF
  240    CONTINUE
  250 CONTINUE
      IF(IDIR.EQ.0) THEN
         CALL pyrobo(nstb,0,STHETAR(2),0.D0,0.D0,0.D0,0.D0)
         CALL pyrobo(nstb,0,0.D0,SPHIR(2),0.D0,0.D0,0.D0)
         CALL pyrobo(nstb,0,0.D0,0.D0,DBETAR(2,1),DBETAR(2,2),
     +   DBETAR(2,3))
      ENDIF
      CHEC = 5.D-4
      P1 =ABS(P(1,1)+P(2,1))
      P2 =ABS(P(1,2)+P(2,2))
      P3 =ABS(P(1,3)+P(2,3))
      P4 =ABS(P(1,4)+P(2,4))
      PC1 = (ABS(pyp(0,1)) - P1)/P4
      PC2 = (ABS(pyp(0,2)) - P2)/P4
      PC3 = (ABS(pyp(0,3)) - P3)/P4
      PC4 = (ABS(pyp(0,4)) - P4)/P4
      IF(DABS(PC1).GT.CHEC.
     +   OR.DABS(PC2).GT.CHEC.
     +   OR.DABS(PC3).GT.CHEC.
     +   OR.DABS(PC4).GT.CHEC) THEN
         write(6,*) ' ppmeps: energy of final particles not correct'
     +   ,chec
         write(6,*) ' PC1 = ',PC1,pC2,pC3,pC4
         write(6,*) ' IPRO = ',IPRO
         write(6,*) 'PYVAR(31),PYVAR(32) ',PYVAR(31),PYVAR(32)
         IF(IDIR.EQ.0) write(6,*) ' M_X = ' ,sqrt(RMX2), ' t = ',T2GKI
         call pylist(1)
         LST(21) = 100
c         stop
      ENDIF
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' in partonic cm system '
         CALL pylist(1)
      ENDIF
c... change z coordinate  to rapgap standard
      do 260 i=1,n
  260 p(i,3) = - p(i,3)
      DO 270 I=1,4
         P(1,I)=PBEAM(1,I)
         P(2,I)=PBEAM(2,I)
  270 CONTINUE
      IF(NCALL.LE.NPRIN) THEN
         write(6,*) ' after boost dbeta in ppmeps'
         CALL pylist(1)
      ENDIF
      IF(NCHECK.EQ.1) THEN
         write(6,*) ' ppmeps final '
         CALL pylist(1)
      ENDIF
      RETURN
      END
