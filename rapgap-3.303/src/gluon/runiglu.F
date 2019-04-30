C*********************************************************************

      SUBROUTINE RUNIGLU(KF,XI,KT2I,TI,XPQ)
      Implicit None
#include "rgparas.inc"
      Integer kf,i,Iccfmset
      REAL XI,KT2I,TI
      REAL XPQ(-6:6)
      Double Precision X,KT2,T,scale
      Double Precision XPQS(-6:6)
      REAL DGRV_NLN
      Real R02, LAM, SIGM0, X0, glu, ALPHAS
      DOUBLE PRECISION DGRV_NL,BLUEML
      Integer IGLU
      Common/gluon/Iglu
      Integer ISEL
      EXTERNAL DGRV_NLN,DGRV_NL,BLUEML
      LOGICAL FIRST
      Data FIRST/.TRUE./
c ISEL = 1 CCFM
c      = 2 derivative of GRV
c      = 3 Bluemlein
c      = 4 KMS
c      = 5 Golec-Biernat Wuesthoff
c      = 6 KMR
      IF(FIRST) THEN
	   If(IGLU.GT.1000) then
	   ISEL = Int(IGLU/1000)
c	   write(6,*) ' isel ',ISEL
	   ELSE
         ISEL = IGLU
	   ENDIF
         FIRST=.FALSE.
         Write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
         Write(6,*) '++++  gluon density selection    ++++++++++++++'
         Write(6,*) '++++  ISEL= ',ISEL,
     +      '                  ++++'
         If(Isel.eq.1) then 
	     if(iglu.lt.1000) then
           Write(6,*) '++++  ISEL=1 --> CCFM   JS2001             ++++'
	     elseif(iglu.gt.1000) then
	        Iccfmset =  MOD(IGLU,1000)
     	        write(6,*) '++++  Iccfmset = ',Iccfmset
	        if(iccfmset.eq.1) then
           Write(6,*) '++++  CCFM  J2003 set 1 (Ipgg0ns-1)        ++++'
	        elseif(iccfmset.eq.2) then
           Write(6,*) '++++  CCFM J2003 set 2 (Ipgg1ns1)          ++++'
	        elseif(iccfmset.eq.3) then
           Write(6,*) '++++  CCFM J2003 set 3 (Ipgg2ns2)          ++++'
	        elseif(iccfmset.eq.10) then 
           Write(6,*) '++++  CCFM set A0                          ++++'
	        elseif(iccfmset.eq.11) then
           Write(6,*) '++++  CCFM set A0+                         ++++'
	        elseif(iccfmset.eq.12) then
           Write(6,*) '++++  CCFM set A0-                         ++++'
	        elseif(iccfmset.eq.13) then
           Write(6,*) '++++  CCFM set A1                          ++++'
	        elseif(iccfmset.eq.20) then
           Write(6,*) '++++  CCFM set B0                          ++++'
	        elseif(iccfmset.eq.21) then
           Write(6,*) '++++  CCFM set B0+                         ++++'
	        elseif(iccfmset.eq.22) then
           Write(6,*) '++++  CCFM set B0-                         ++++'
	        elseif(iccfmset.eq.23) then
           Write(6,*) '++++  CCFM set B1                          ++++'
	        elseif(iccfmset.eq.101) then
           Write(6,*) '++++  CCFM set C                           ++++'
	        elseif(iccfmset.eq.111) then
           Write(6,*) '++++  CCFM test set                        ++++'
	        else
	     Write(6,*) '++++ no valid set specified:',iglu,'++++'
	        endif
	     endif
	     elseif(isel.eq.2) then
         Write(6,*) '++++  ISEL=2 --> derivative of coll. gluon ++++'
	     elseif(isel.eq.3) then
         Write(6,*) '++++  ISEL=3 --> Bluemlein                 ++++'
	     elseif(isel.eq.4) then
         Write(6,*) '++++  ISEL=4 --> Kwiecinski,Martin,Stasto  ++++'
	     elseif(isel.eq.5) then
         Write(6,*) '++++  ISEL=5 --> Wuesthoff Golec-Biernat   ++++'
	     elseif(isel.eq.6) then
         Write(6,*) '++++  ISEL=6 --> Kimber,Martin,Ryskin      ++++'
	     elseif(isel.eq.7) then
         Write(6,*) '++++  ISEL=7 --> Ryskin,Shabelski          ++++'
	   endif
         Write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
      ENDIF

      X = XI
      KT2 = KT2I
c for new p_max
      T=max(0.0001,TI)
c     T=0.
c unintegrated gluon for semihard approach
      DO 10 I=-6,6
   10 XPQ(I)=0.0
      SCALE = KT2
      If(ISEL.EQ.1) THEN
	   if(iabs(kf).eq.2212) then 
            call ccfmstfu(kf,x,scale,t,xpqs)
c            write(6,*) xpq(0)
	   else
            write(6,*) ' cauniglu: no valid kf = ',kf	      
	   endif 

c here comes the test stuff
      ELSEIF(ISEL.EQ.2) THEN
         SCALE = MAX(0.35,KT2)
c analytical derivative of grv (fast)
         XPQs(0) = dgrv_nln(X,SCALE)
c numerical derivative of gluon (slow)
c             XPQ(0) = d_xgx(X,SCALE)
         IF(XPQs(0).LE.0.) XPQs(0)=0.
      ELSEIF(ISEL.EQ.3) THEN
c select Bluemlein parametrisation
         XPQs(0) = SNGL(BLUEML(X,SCALE,T))
      ELSEIF(ISEL.EQ.4) THEN
         SCALE = KT2
         if(kt2.lt.1.) then
            scale = 1.01
            call ktstfu(kf,x,scale,xpqs)
         else
            if(kt2.gt.100.) scale = 100.
            call ktstfu(kf,x,scale,xpqs)
         endif
      ELSEIF(ISEL.EQ.5) THEN
         SCALE = KT2
cccc this is from SATRAP
c ori      SIGM0 = 23.03*2.568
c          mb -> 1/GeV^2 = 1/0.389379 = 2.568
c         SIGM0 = 23.03d0*2.568d0 ! w/o charm
*! with charm
         SIGM0 = 29.12d0*2.568d0
c ori      LAM  = 0.288  ! w/o charm
*! with charm
         LAM = 0.277d0
c ori      X0 =  3.04/10000. ! w/o charm
*! with charm
         X0 = 0.41d0/10000.d0
         ALPHAS = 0.2d0
         R02 = (dble(X)/X0)**LAM
         glu = 3.d0*SIGM0*R02*SCALE*EXP(-R02*SCALE)
	   if(ipro.eq.20) then
           glu = 3.d0/(2.25*2.25)*SIGM0*R02*dble(SCALE)*
     &     EXP(-R02*dble(SCALE)/2.25)
	   endif
c pi**2 = 9.8696
         XPQS(0) = glu/4.d0/9.8696d0/ALPHAS
      ELSEIF(ISEL.EQ.6) THEN
	   if(kf.eq.22.or.iabs(kf).eq.11) then
	   write(6,*) ' cauniglu : '
	   write(6,*) ' requested set',isel,' not available for kf = ',kf
         Stop
	   endif 
         call ccfmstfu(kf,x,scale,t,xpqs)
      ELSEIF(ISEL.EQ.7) THEN
	   if(kf.eq.22.or.iabs(kf).eq.11) then
	   write(6,*) ' cauniglu : '
	   write(6,*) ' requested set',isel,' not available for kf = ',kf
         Stop
	   endif 
         call ryssha(kf,x,scale,xpqs)
      ELSE
         write(6,*) ' ISEL = ',ISEL,' not implemented '
         STOP
      ENDIF
      DO I=-6,6
         XPQ(I)=XPQS(I)
      ENDDO
      RETURN
      END
