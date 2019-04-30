      SUBROUTINE ELQQFNL(WT1)
C  calculates the matrixelement 
C for light quarks without masses using the subtraction scheme
c changed for including Q2 of photon
C
      IMPLICIT NONE
#include "rgpart.inc"
#include "rglujets.inc"
#include "rgpara1.inc"
#include "rgpara.inc"
#include "rgluco.inc"
#include "rgrapgki.inc"
#include "rgdisdif.inc"
#include "rgnlo.inc"
      Double Precision WMAX
	Integer IMIX
      COMMON /OALPHAS/ WMAX,IMIX
      Double Precision P2E(4),Q1E(4),Q2E(4)
      Double Precision ALPHA_S,ALPHAS,ALPH_EM,SH,TH,UH,QF2
      Double Precision EPSILON,YX,WT1,DOT
      Double Precision p_gqq, Z, Z1,Z2,q 
      Double Precision X1,X1A
      double precision antiux1
      double precision antiqu
      double precision ux1
      double precision ux,antiux,xbj
      double precision ctnlo,const2
      double precision up,hq,restbgf
      double precision antiup,antihq
      double precision qu,xold
      double precision scal,cost,x3,cosfac
	Double Precision ME,fgam
	Integer I
      REAL SNGL
      REAL XPQ(-6:6)
	External ctnlo
      Logical NLO,INi_NLO
      Integer Nscheme, Nfmx, ifil
      Common /NLOpara/ NLO,Ini_nlo, NSCHEME, NfMx, ifil
	Logical ini_old

      ini_old=ini_nlo
	ini_nlo=.true.
c we set ini-nlo=true to obtain the original msbar pdfs  
c via the call to rastfu/rgstfu	
C------------------RESET-------------------
      ux=0.0d0
      antiux=0.0d0
      ux1=0.0d0
      antiux1=0.0d0
      wt1=0.0d0
	hq=0d0
	antihq=0d0
	restbgf=0d0
      DO I = -6, 6
         XPQ(I) = 0.0
      enddo
	ME = PYMASS(K(1,1))
C------------------------------------------
cccc      if (kpa.ge.0) then  ! sabine ?????
      IF(INTER.EQ.2) THEN
         write(6,*) ' interaction INTER = 2 not implemented '
         STOP
      ENDIF
      DO 10  I=1,4
         P2E(I) = P(NIA2,I)
	   IF(KPA.GT.0) then
           Q1E(I) = P(NF1,I)
           Q2E(I) = P(NF2,I)
	     else
           Q1E(I) = P(NF1,I)
           Q2E(I) = P(NF2,I)
	   endif
   10 CONTINUE
   20 CONTINUE
C     MARTIX ELEMENT
      if (q2q.ne.q2) then
         print*,'warning: the subtraction scheme can be only'
         print*,'used at the scale mu**2=Q**2'
         print*,'q2,q2q in ellleqql',q2,q2q
         stop
      endif
      ALPHA_S=ALPHAS(DSQRT(Q2))
c      alphas2=ALPHAS(DSQRT(Q2q))

      SH =  2.D0 * DOT(Q1E,Q2E)
      TH = - 2.D0 * DOT(P2E,Q2E)
      UH = - 2.D0 * DOT(P2E,Q1E)
c sum of charges for u d s (2/3)**2 + (1/3)**2 + (1/3)**2 = 6/9 = 2/3
      
c       QF2=2.D0/3.D0
      IF(IMIX.EQ.1) QF2 = DFLOAT(PYCHGE(KPA))**2/9.D0
c      print*,'qf2',qf2,imix
c epsilon
      YX = DBLE(YY)
      EPSILON = (1.D0 - YX)/(1.D0 - YX + YX**2 /2)
c    epsilon end 
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = PYALEM(SNGL(Q2))
c      CONST=16.D0*PI*PI*ALPHA_S*ALPH_EM*QF2
C     NOW WE HAVE TO DEFINE THE SHOWER CONTRIBUTION 
c     FIRST: KINEMATICS
      XBJ=Q2/(SSS*YX)
      IF(IDir.EQ.1) THEN
      ELSE
         XOLD=XBJ
         XBJ=XOLD/XFGKI
         IF(XBJ.LT.1.and.XBJ.gt.0.) THEN
c            write(6,*) ' pqcd beta found,xfgki',xold,xfgki,xbj
         ELSE
            goto 123
         ENDIF
      endif
c         print*,'idir,xbj,xold', idir,xbj,xold/xfgki
c     print*,'xbj',xbj
      q=dsqrt(Q2)
      z=Q2/(SH+Q2)
      x3=xbj/z
c         print*,'x3,xpr',x3,xpr
      cost=((-2.d0*z*uh/q2) -1.0d0)
	if(abs(cost).gt.1d0) then
	   write(6,*) ' cost > 1 ',cost
	   write(6,*) ' z,uh,q2 ',z,uh,q2 
	endif
      Z1=Z+TH*Z/Q2
      Z2=Z+UH*Z/Q2
      x1=xbj*(1.d0+th/Q2)   
      x1a=xbj*(1.d0+uh/Q2) 
	if(abs(cost).gt.0.999) goto 123
c     the right scale is always q2
      scal=q2
c     idir =1: normal DIS xbj
c     idir =0: Diffractive case beta=xbj/xpom
      IF(IDIR.EQ.1) THEN
         CALL rgstfu(2212,SNGL(XBJ),SNGL(SCAL),XPQ)
         ELSE
         CALL RASTFU(KINT(2,2),SNGL(XBJ),SNGL(SCAL),XPQ)
      ENDIF
      ux=dble(XPQ( kpa))/xbj
      antiux=dble(XPQ(-kpa))/xbj
C     QUARK CONTRIBUTION
      if(ctnlo(q2,-th).gt.0d0) then
	   cosfac = 1d0 - cost
         if(th.eq.0.0d0.or.cost.eq.1.d0) then 
            hq=0.0d0
            print*,'t,cost,q2 hardquark',th,cost,q2,hq
         else 
	      if(x1.lt.0) then
		  write(6,*) ' elqqfnl x1 < 0 ',x1,th,uh,q2,cost
		endif
            IF(IDIR.EQ.1) THEN
               CALL rgstfu(2212,SNGL(X1),SNGL(SCAL),XPQ)
            ELSE
               CALL RASTFU(KINT(2,2),SNGL(x1),SNGL(SCAL),XPQ)
            ENDIF
            ux1=dble(XPQ(kpa))/x1
            if(ux1.eq.0.0d0) then !only antiux1 becomes 0, rest is security
               qu=0.d0
            else
               qu=ux/ux1
c		   qu = 1d0
            endif 
ccccold            up=(p_gqq(Z) -ctnlo(q2,-th)*p_gqq(Z1) *qu)/cosfac
ccc		up = 0d0
		up = p_gqq(Z)*(1d0 - ctnlo(q2,-th)*qu) 
		if(abs(up).lt.0.0000001) then
		   up = 0.d0 
		   else 
		   up = up/cosfac
		endif
		up = up - (1d0 - 2d0*z + cosfac/2d0)*ctnlo(q2,-th)*qu/2d0
            hq=up
c        print*,'th,cost,hardquark,up',th,cost,hq,up
c     print*,'z,z1,up,down,strange',z,z1,up
         endif 
	else
c	   write(6,*) ' ctnlo(q2,-th) = 0 ,q2,th' ,ctnlo(q2,-th),q2,th
	   cosfac = 1d0 - cost
c	   write(6,*) ' 1 - cost ',cosfac
	   hq=p_gqq(Z)/cosfac
	endif
*************
      if(ctnlo(q2,-uh).gt.0d0) then
	   cosfac = 1 + cost
         if(uh.eq.0.0d0.or.cost.eq.-1.d0) then
            antihq=0.0d0
            print*,'u,q2,cost, antihardquark',uh,q2,cost,antihq
         else
	      if(x1a.lt.0) then
		  write(6,*) ' elqqfnl x1a < 0 ',x1a,th,uh,q2,cost
		endif
            IF(IDIR.EQ.1) THEN
               CALL rgstfu(2212,SNGL(X1a),SNGL(SCAL),XPQ)
            ELSE
               CALL RASTFU(KINT(2,2),SNGL(x1a),SNGL(SCAL),XPQ)
            ENDIF
            antiux1=dble(XPQ(-kpa))/x1a
            if(antiux1.eq.0.0d0) then !only antiux1 becomes 0, rest is security
               antiqu=0.d0
            else
               antiqu=antiux/antiux1
c               antiqu=1d0
            endif 
ccc            antiup=(p_gqq(Z) -ctnlo(q2,-uh)*p_gqq(Z2)* antiqu)/cosfac
ccc            antiup=0d0
		antiup = p_gqq(Z)*(1d0 - ctnlo(q2,-uh)*antiqu)
		if(abs(antiup).lt.0.0000001) then
		   antiup = 0.d0 
		   else 
		   antiup = antiup/cosfac
		endif
		antiup = antiup 
     &		-(1d0-2d0*z+cosfac/2d0)*ctnlo(q2,-uh)*antiqu/2d0
            antihq=antiup
c        print*,'uh,cost,antihardquark,antiup',
c     &     uh,cost,antihq,antiup
         endif
	else
c	   write(6,*) ' ctnlo(q2,-uh) = 0 ,q2,uh' ,ctnlo(q2,-uh),q2,uh
	   cosfac = 1d0 + cost
c	   write(6,*) ' 1 + cost ',cosfac
	   antihq=p_gqq(Z)/cosfac
	endif

      restbgf= 3.d0*z*(1d0-z)-1.d0/2.d0
c the following factor corresponds to the factors in eq(3.21) of J.Collins
c the additional factor 8 comes from the phase space in PHASE.F which is:
C dLips = 1/8 dcos dphi	
      CONST2 = 8.d0*alpha_s*alph_em**2*qf2*xbj/x3/q2/pi
	fgam = ((1.d0-yx+yx*yx/2d0)/q2/yx - me*me*yx/q2/q2)
	const2=const2*fgam
	wt1 = hq+antihq+restbgf
	if(wt1.lt.0d0) then
	   if(abs(wt1).lt.0.0001) wt1=0d0 
	endif
      wt1=const2*wt1
	ntotnlo=ntotnlo + 1
	if(wt1.ne.wt1) then
	   write(6,*) ' elqqfnl problem wt1 NAN ',wt1
	   wt1=0d0
	   elseif(wt1.lt.0d0) then
	   write(6,*) ' elqqfnl problem wt1 < 0 ',wt1
	   write(6,*) ' hq+antihq+restbgf ',hq+antihq+restbgf
	   write(6,*) ' hq,antihq,restbgf,cosfac',hq,antihq,restbgf,
     &    cosfac
	   write(6,*) ' up, antiup',up, antiup
         write(6,*) ' p_gqq(Z)',p_gqq(Z)
	   write(6,*) ' p_gqq(Z1) *qu',p_gqq(Z1)*qu
	   write(6,*) ' p_gqq(Z2) *qu',p_gqq(Z2)*antiqu
	   negnlo=negnlo + 1
	   wt1=0d0
	endif
123   Continue
	ini_nlo=ini_old
c	write(6,*) ' wt1 , cosfac ',wt1,cost    
      RETURN
      END
      
      
    

      
