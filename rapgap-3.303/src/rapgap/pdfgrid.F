	Subroutine pdfgrid
	Implicit None
#include "rgnlo.inc"
#include "rgdiffr.inc"
#include "rgpart.inc" 
#include "rgluco.inc"
#include "rgpara1.inc"
#include "rgparas.inc"
#include "rgpypars.inc"
#include "rglujets.inc"
#include "rgnlogrid.inc"
#include "rgdisdif.inc"
#include "rgraher.inc"
#include "rgnlodiff.inc"
#include "rgrapgki.inc"
	Double Precision XXI(NPDFX),Q2I(NPDFQ2)
      double precision xpom
	Integer ID
      COMMON /PQPMID/xpom,ID
		
	double Precision xpq(-6:6)
	Double Precision dq,dx,x,q,xmif,q2mif,q2maf
	Integer mstp51,mstp51f,mstp52,mstp52f,icutself,idirf
	integer I,J,II,ik
	Integer iostat,irr,err,ifil
	Real xpqdum(-6:6)
	
	Double Precision QiniN, QmaxN, XminN, XcrN
      Common / NLOXQrange / QiniN, QmaxN, XminN, XcrN
      Logical Ini_NLO, NLO
	Integer Nscheme, NfMx
      Common /NLOpara/ NLO, INI_NLO,NSCHEME, NfMx, ifil
	Double Precision NLOpdf
	External NLOpdf
	logical ex,lprint
c	Data lprint/.false./
	Data lprint/.true./
	
c     parameters for NLOset
      NFMX=5
	xcrN = 1.5d0
      QiniN=2d0 
	QmaxN=1.0d+4 
	xminN=1.0d-4
      Ini_NLO = .True.
	IHF = 0
	ISEMIH = 0
	err = 0
	iostat = 0
	
	ID = 1 
      IF(IDIR.EQ.1) THEN
         ID = 1
      ELSEIF(IDIR.EQ.0.AND.(NPOM.EQ.20.OR.NPOM.EQ.21)) THEN
         ID = 3
      ELSE
         ID = 2
      ENDIF
C select leading order or higher order process
      IF(IDISDIF.GE.3.AND.IDIR.EQ.1) THEN
         write(6,*) ' this option is  not impelmented --> STOP '
         STOP
      ENDIF
      IF(IDISDIF.GE.1.AND.IDIR.EQ.1) THEN
         IHERPYS = 0
         ID = 1
      ENDIF
      IF(ID.EQ.1) THEN
            CALL rgstfu(2212,0.01,10.,XPQDUM)
		Diff_flux =1.
            Q2MIN = 0.50D0
c            Q2MAX = DMIN1(QMA,SSS)
            Q2MAX = SSS
c            XMIN = SNGL(Q2MIN/SSS)
            XMIN = 1.D-06  ! for non-diff
         ELSE
            T2GKI = -0.1
            XFGKI = 0.01
            CALL RASTFU(KINT(2,2),0.01,10.,XPQDUM)
            Q2MIN = 4.00D0
            Q2MIN = 0.50D0
c            Q2MAX = DMIN1(QMA,SSS)
            Q2MAX = SSS
c            XMIN = SNGL(Q2MIN/SSS)
            XMIN = 1.D-05  ! for non-diff
         ENDIF
	
	
      DQ = (LOG(Q2MAX) - LOG(Q2MIN))/DBLE(NPDFQ2-1)
      DX = (LOG10(0.999D0) - LOG10(XMIN))/DBLE(NPDFX-1)

      QiniN = dsqrt(q2min)
	QmaxN = dsqrt(q2max)
	XminN = XMIN-1.d-10

	MSTP51=MSTP(51)
	MSTP52=MSTP(52)

	write(6,*) 
     &' Calculate conversion from MSbar to BS parton densities ' 
      write(6,*) ' with  xmin = ',xmin
	write(6,*) ' with q2min = ',q2min,' q2max = ',q2max  
             
      inquire(FILE='bspdf.dat',EXIST=ex)
      if(ex) then
         open(32,FILE='bspdf.dat', FORM='unformatted',STATUS= 'OLD',
     +   IOSTAT=IRR,ERR=220)
         read(32) mstp51f,mstp52f,icutself,idirf
         read(32) xmif,q2mif,q2maf
         if(xmif.eq.xmin.and.q2mif.eq.q2min.and.
     &   q2maf.eq.q2max.and.mstp51f.eq.mstp51.and.mstp52f.eq.mstp52.and.
     &   icutself.eq.icutsel.and.idirf.eq.idir) then
            write(6,*) ' bspdf.dat found with consistent parameters '
		ifil = 1
            else
		close(32)
            write(6,*) ' inconsistent parameters on file on unit '
     +      //'32 bspdf.dat '
            open(32,FILE='bspdf.dat', FORM='unformatted',STATUS=
     +      'OLD',IOSTAT=IRR,ERR=220)
		ifil = 0
            write(6,*) ' file found, overwriting old file bspdf.dat'
            write(32) mstp51,mstp52,icutsel,idir 
            write(32) xmin,q2min,q2max
         endif
         else
         open(32,FILE='bspdf.dat', FORM='unformatted',STATUS=
     +   'NEW',IOSTAT=IRR,ERR=220)
         write(6,*) ' no file found on unit 32 '
         write(6,*) ' cold start necessary '
         ifil = 0
         write(32) mstp51,mstp52,icutsel,idir
         write(32) xmin,q2min,q2max
	endif

	Call NLOset 
	
      if(LPRINT) write(6,9999)
      DO I=1,NPDFX 
         DO J=1,NPDFQ2-1
         XXI(I) = SNGL(10**(LOG10(xmin) + DX*DFLOAT(I-1)))
         IF(XXI(I).GE.1.D0) XXI(I) = 0.9999
         Q2I(J) = SNGL(EXP(LOG(Q2MIN) + DQ*DFLOAT(J-1)))
         XPDF(I) = DBLE(XXI(I))
         Q2PDF(J) = DBLE(Q2I(J))
         x = xxi(i)
         q = dsqrt(q2i(j))
c	   write(6,*) ' pdfgrid x ',x,q**2,sss
c for the light flavors	   
         DO II = -6,6
            XPQ(II)=sngl(NLOpdf(ii, x, q))
         enddo
	   IF(LPRINT) then
	   write(6,1000) x,q,(xpq(ik),ik=-3,3)
	   if(ID.EQ.1) THEN
            CALL rgstfu(2212,sngl(xpdf(i)),sngl(q2pdf(j)),XPQDUM)
	   else
            CALL RASTFU(KINT(2,2),sngl(xpdf(i)),sngl(q2pdf(j)),XPQDUM)
         ENDIF
	   write(6,1000) x,q,(xpqdum(ik)/diff_flux,ik=-3,3)
	   endif
      enddo
	enddo
	IGRICALC = 1
      Ini_NLO = .False.
	close(32)
	Return
9999  FORMAT('    X     |   Q     |   sbar  |',
     +'  ubar  |  dbar  |  gluon |   d    |   u    |    s   |')
 1000 FORMAT(' ',E8.3,' |',E8.3,' | ',7(E8.3,'|'))
220   write(6,*) ' error opening bspdf.dat file '
	stop
	End
      
      double precision function intquark(x,q2)
      implicit none
      double precision x,q2,xi,q2i
      real sumquark,quark,eps
      common/insidepart/xi,q2i
      external fl2,fu2,quark
	
      xi=x
      q2i=q2
	
      eps=0.01
      call gadap2 (0.0,1.0,fl2,fu2,quark,eps,sumquark)
      if(sumquark.le.0) then
c	   write(6,*) ' intquark: quark density < 0: ',sumquark,
c     &	   ' at ',x,q2
         endif
      intquark =dble(sumquark)
c     print*,'sumquark,x,q2,intoquark,',x,q2,sumquark,intquark
      return
      end
	
      function quark(R,rcost)
      implicit none
#include "rgdisdif.inc"
#include "rgpart.inc"
#include "rgnlodiff.inc"
      double precision x,q2,p_gqq,xi,q2i
      double precision x1,z,x3,cost
      double precision mu,mu2,bracket,bracketn
      double precision fa,fa1,th,uh,virt,cosfac,co1,glu,dzglu
	Double Precision ctnlo
      real rcost,rx3,quark,weight,R,xmin,xmax
      external p_gqq
      common/insidepart/xi,q2i
	Integer ifl
      common/flavour/ifl
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      REAL SNGL
      REAL XPQ(-6:6)
	External ctnlo
      x=xi
      q2=q2i
	xmin=x
	xmax=1.
      rx3=(xmin*(xmax/xmin)**R)
      weight=rx3*log(xmax/xmin)
c      print*,'x in quark ',x
      cost=dble(rcost)
      x3=dble(rx3)
      mu=dsqrt(q2)
	mu2=mu**2
      z=x/x3
c	write(6,*) ' before quark stfu call '
      IF(ID.EQ.1) THEN
            CALL rgstfu(2212,RX3,SNGL(MU2),XPQ)
         ELSE
            CALL RASTFU(KINT(2,2),RX3,SNGL(MU2),XPQ)
         ENDIF
c	write(6,*) ' aftr quark stfu call '
      if(xpq(0).ne.xpq(0)) then
	  write(6,*) ' xpq(0) NAN ... ',xpq(0),rx3,mu2
	  endif
      glu=DBLE(XPQ(0))/diff_flux/x3
c	write(6,*) ' glu   ',xpq(0),rx3,diff_flux
	if(glu.gt.0) then
      dzglu=z*glu 
      co1=(p_gqq(z)*dlog((Q2*(1.d0-z))/(mu**2*z)) +z*(1.d0-z))/2.d0  
C   divided by two to compensate
C   the contribution from the cost--integration
      th = -Q2*(1.d0-cost)*x3/(2.0d0*x)
      uh = -Q2*(1.d0+cost)*x3/(2.0d0*x)
c select proper virtuality, t or u, accordning to quark or antiquark
      if(ifl.gt.0) then
	   virt = abs(th)
	   cosfac = 1d0 - cost
         x1=x*(1.0d0+th/Q2) 
	else
	   virt = abs(uh)
	   cosfac = 1d0 + cost
         x1=x*(1.0d0+uh/Q2) 
	endif
	IF(ctnlo(q2,virt).gt.0d0.and.x1.le.0) then
	   write(6,*) ' problem in quark: ctnlo = ',ctnlo(q2,virt)
	   write(6,*) ' virt = ',virt,' q2 = ',q2,' x1 = ',x1
	endif
C     x1test=x-0.5d0*x3*(1.d0-cost)
c     print*,'x,x1,x3',x,x1,x1test,x3
	IF(ctnlo(q2,virt).gt.0d0) then
	   IF(ID.EQ.1) THEN
            CALL rgstfu(2212,SNGL(X1),SNGL(MU2),XPQ)
         ELSE
            CALL RASTFU(KINT(2,2),SNGL(x1),SNGL(MU2),XPQ)
         ENDIF
         fa1= dble(XPQ(ifl))/diff_flux/x1
c     define the kinematics as in 3.5 Collins in order to avoid the commonblock!   
         if(FA1.EQ.0.0d0) then
            bracketn=0.d0
         else 
            IF(ID.EQ.1) THEN
               CALL rgstfu(2212,SNGL(X),SNGL(MU2),XPQ)
            ELSE
               CALL RASTFU(KINT(2,2),SNGL(x),SNGL(MU2),XPQ)
            ENDIF
            fa=dble(XPQ(ifl))/diff_flux/x
C                                ! msbar at mu**2          
c            z1=z-0.5d0*cosfac
c            bracket=p_gqq(z)-ctnlo(q2,virt)*fa/fa1*p_gqq(z1)
c            bracketn=bracket/cosfac

		bracket = p_gqq(Z)*(1d0 - ctnlo(q2,virt)*fa/fa1) 
		if(abs(bracket).lt.0.0000001) then
		   bracket = 0.d0 
		   else 
		   bracket = bracket/cosfac
		endif
c		write(6,*) ' bracket 1 ',bracket
		bracket = bracket - 
     &      (1d0 - 2d0*z + cosfac/2d0)*ctnlo(q2,virt)*fa/fa1/2d0
c		write(6,*) ' bracket 2 ',bracket
            bracketn=bracket
		
            if(FA1.EQ.0.) then 
               print*,'ifl,x,bracketn in quark ',ifl,x,bracketn
            endif
            endif
	   else
            bracket=p_gqq(z)
            bracketn=bracket/cosfac
c	   write(6,*) ' no sub', bracketn,bracket,cosfac
	   endif
         quark = real(dzglu*(co1-bracketn)*weight)
c	   if(quark.le.0) then 
c	     write(6,*) ' quark < 0: quark = ',quark,' z= ',z
c	     write(6,*) ' quark < 0: co1 ',co1,' bracketn ',bracketn
c	     write(6,*) ' quark < 0: th,uh,q2  ',th,uh,q2
c	     write(6,*) ' quark < 0: glu  ',glu,x3,z,cost
c           quark=0
c	   endif
	   if(quark.ne.quark) then 
	      write(6,*) ' quark NAN: quark = ',quark,dzglu,co1,bracketn
	      write(6,*) ' quark NAN: th,uh,q2  ',th,uh,q2
	      write(6,*) ' quark NAN: glu  ',glu,x3,z,cost
         endif
	else
	   quark = 0d0
	endif
      return
      end

*******************
      
      function fu2(rcost)
      implicit none
      real fu2,rcost
      
c      fu2=1.0
      fu2=0.9999
      return
      end 

******************
      function fl2(Rcost)
      implicit none
      real fl2,rcost
c      fl2=-1.
      fl2=-0.9999
      return 
      end
*
C==================================================================
*****************************************************************************

C March 23, 2002, H. Jung

C july 29, 2001, Xiaomin
C NLOtrans, given x, q, and Iprtn ( parton flavor), this function
C transform msbar pdf to NLOpdf, right now only BS scheme. 
C Modified from Sabine Schilling's code bs.f
C
*****************************************************************************
      Function NLOtrans(Iprtn, x, q)
      implicit none
	Integer Iprtn 
	Double precision NLOtrans
      double precision alpha_s,q,x
      double precision xfa,xfabs,p_gqq,q2i
      double precision intquark,xi
      double precision alphas
#include "rgdisdif.inc"
#include "rgparam.inc"
#include "rgpara1.inc"
#include "rgpart.inc"
#include "rgnlodiff.inc"
      external p_gqq
      external intquark
	Integer ifl
      common/flavour/ifl
      common/test/alpha_s
      common/insidepart/xi,q2i
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      double precision beta,scale2
      common /spom/ beta,scale2

      REAL SNGL
      REAL XPQ(-6:6)
      
      q2=q**2
      q2i=q2
c     beta=x
      xi=x
c     print*,'x in NLOtrans',x
	NLOtrans=0d0
c	write(6,*) ' start: NLO trans ',xfabs,x,q2
      pi = 4.d0*atan(1.d0)
      alpha_s=alphas(q)
c	write(6,*) ' before Nlotrans stfu call '

      IF(ID.EQ.1) THEN
         CALL rgstfu(2212,SNGL(X),SNGL(q2),XPQ)
         ELSE
         CALL RASTFU(KINT(2,2),SNGL(x),SNGL(q2),XPQ)
      ENDIF
c	write(6,*) ' after Nlotrans stfu call '
      ifl = Iprtn
	xfabs = 0d0
	if(iabs(ifl).le.3) then
c	  write(6,*) ' NLO trans in loop '
        xfa=dble(XPq(ifl))/diff_flux  
	  If(ifl.ne.0) then
c     fabs is parton densities in BS scheme at x: 
           xfabs=xfa + alpha_s/(2.d0*pi)*intquark(x,Q2)
c           xfabs=xfa  
	     else
	     xfabs=xfa
	  endif
        if (xfabs.le.0.d0) then 
           xfabs=0.0d0
c          print*,' the parton density is set to 0, as negative value',
c    &     xfabs,xfa,alpha_s/(2*pi)*intquark(x,Q2),dxbs(a)
        endif 
      endif
c	write(6,*) ' NLO trans ',xfabs,x,q2
      NLOtrans = xfabs
      return
      end

      Function NLOtmp(Iprtn, x, q)
      implicit none
	Integer Iprtn 
	Double Precision NLOtmp
      double precision alpha_s,q,x,mu,pymass
      double precision xfa,xfabs,p_gqq,q2i
      double precision inttmp,intQQ,xi,bgf
      double precision alphas
#include "rgdisdif.inc"
#include "rgparam.inc"
#include "rgpara1.inc"
#include "rgpart.inc"
#include "rgnlodiff.inc"
      external p_gqq
      external inttmp,intQQ,pymass
	Integer ifl
      common/flavour/ifl
      common/test/alpha_s
      common/insidepart/xi,q2i
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      double precision beta,scale2
      common /spom/ beta,scale2

      REAL SNGL
      REAL XPQ(-6:6)
      
      q2=q**2
	
      q2i=q2
	
c     beta=x
      xi=x
c     print*,'x in NLOtmp',x
	NLOtmp=0d0
c	write(6,*) ' start: NLO tmp ',xfabs,x,q2
      pi = 4.d0*atan(1.d0)
      alpha_s=alphas(q)
      IF(ID.EQ.1) THEN
         CALL rgstfu(2212,SNGL(X),SNGL(q2),XPQ)
         ELSE
         CALL RASTFU(KINT(2,2),SNGL(x),SNGL(q2),XPQ)
      ENDIF
      ifl = Iprtn
	If(Iabs(ifl).le.3) then
        xfa=dble(xpq(ifl))/diff_flux  
	  If(ifl.ne.0) then
c     fabs is parton densities in BS scheme at x: 
           bgf=inttmp(x,Q2)
           xfabs=xfa + alpha_s/(2.d0*pi)*bgf
ccc           xfabs=xfa  
	     else
	     xfabs=xfa
	  endif
        if (xfabs.le.0.d0) then 
	     write(6,*) ' problem in calcualting NLO BGF contribution'
	     write(6,*) ' xfa = ',xfa,' BGF ',bgf,' ifl =',ifl
	     write(6,*) ' x = ',x,' mu = ',q2
c           xfabs=0.0d0
c          print*,' the parton density is set to 0, as negative value',
c    &     xfabs,xfa,alpha_s/(2*pi)*intquark(x,Q2),dxbs(a)
        endif 
      else
c here calculate heavy quark contribution
       mu=pymass(ifl)
       alpha_s=alphas(2.d0*mu)
       xfabs=          alpha_s/(2.d0*pi)*intQQ(x,Q2)
	endif
c	write(6,*) ' NLO tmp ',xfabs,x,q2
      NLOtmp = xfabs
      return
      end
      double precision function inttmp(x,q2)
      implicit none
      double precision x,q2,xi,q2i
      real sumquark,deltams,eps
      common/insidepart/xi,q2i
	external deltams
      xi=x
      q2i=q2
	
      eps=0.01
c      print*,' before gadap ',x,q2 
      call gadap(0.0,1.0,deltams,eps,sumquark)
      if(sumquark.le.0) then
c	   write(6,*) ' intquark: quark density < 0: ',sumquark,
c     &	   ' at ',x,q2
         endif
      inttmp =dble(sumquark)
c      print*,' inttmp: x,q2',x,q2,sumquark
      return
      end
      double precision function intQQ(x,q2)
      implicit none
      double precision x,q2,xi,q2i
      real sumquark,bgfcc,eps
      common/insidepart/xi,q2i
	external bgfcc
      xi=x
      q2i=q2
	
      eps=0.01
c      print*,' before gadap ',x,q2 
      call gadap(0.0,1.0,bgfcc,eps,sumquark)
      if(sumquark.le.0) then
c	   write(6,*) ' intQQ: quark density < 0: ',sumquark,
c     &	   ' at ',x,q2
         endif
      intQQ =dble(sumquark)
c      print*,' intQQ: x,q2',x,q2,sumquark
      return
      end
	
      function deltams(R)
      implicit none
#include "rgdisdif.inc"
#include "rgpart.inc"
#include "rgnlodiff.inc"
      double precision x,q2,p_gqq,xi,q2i
      double precision z,x3
      double precision mu2
      double precision co1,glu,dzglu
      real rx3,deltams,weight,R,xmin,xmax
      external p_gqq
      common/insidepart/xi,q2i
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      REAL SNGL
      REAL XPQ(-6:6)
      x=xi
      q2=q2i
	xmin=x
	xmax=1.
      rx3=(xmin*(xmax/xmin)**R)
      weight=rx3*log(xmax/xmin)
	deltams = 0.
      x3=dble(rx3)
	mu2=q2
      z=x/x3
c	write(6,*) ' before deltams stfu call '
      IF(ID.EQ.1) THEN
         CALL rgstfu(2212,SNGL(X3),SNGL(MU2),XPQ)
      ELSE
         CALL RASTFU(KINT(2,2),SNGL(X3),SNGL(MU2),XPQ)
      ENDIF
c	write(6,*) ' after deltams stfu call '
      glu=DBLE(XPQ(0))/diff_flux/x3
	if(glu.gt.0) then
         dzglu=z*glu 
         co1=p_gqq(z)*dlog((Q2*(1.d0-z))/(mu2*z)) +
     &	   4.0d0*z*(1.d0-z)-0.5d0
         deltams=dzglu*co1
	else
	   deltams=0d0
	endif
c	if(deltams.le.0) deltams=0
c	write(6,*) ' deltams: z ',deltams,z
	deltams = deltams*weight
c	write(6,*) ' deltams: x3,x,z ',x3,x,z
c	write(6,*) ' deltams: glu,dzglu ',glu,dzglu
c	write(6,*) ' deltams: co1,deltams,weight ',co1,deltams,weight
      return
      end
      function bgfcc(R)
      implicit none
#include "rgdisdif.inc"
#include "rgpart.inc"
#include "rgnlodiff.inc"
      double precision x,q2,xi,q2i
      double precision z,x3
      double precision mu2
      double precision glu,dzglu
	Double Precision am2,PYMASS,beta,zlog,xksi 
      real rx3,bgfcc,weight,R,xmin,xmax
      common/insidepart/xi,q2i
	Integer ifl
      common/flavour/ifl
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      REAL SNGL
      REAL XPQ(-6:6)
      x=xi
      q2=q2i
	am2=PYMASS(ifl)**2
c	write(6,*) ' bgfcc am2 = ',am2,ifl
	mu2=4d0*am2
	xksi=am2/q2
	xmin=x*(1d0+4d0*xksi)
	xmax=1.
	bgfcc=0d0
	if(xmin.ge.xmax) return
      rx3=(xmin*(xmax/xmin)**R)
      weight=rx3*log(xmax/xmin)
	bgfcc = 0.
      x3=dble(rx3)
      z=x/x3
	beta=1d0-xksi*z/(1d0-z)
	beta=dsqrt(beta)
	if(beta.ge.1) return
c	write(6,*) ' before bgfcc stfu call '
      IF(ID.EQ.1) THEN
         CALL rgstfu(2212,SNGL(X3),SNGL(MU2),XPQ)
      ELSE
         CALL RASTFU(KINT(2,2),SNGL(X3),SNGL(MU2),XPQ)
      ENDIF
c	write(6,*) ' after bgfcc stfu call ',xpq(0),mu2,xmax,xmin
      glu=DBLE(XPQ(0))/diff_flux/x3
	if(glu.gt.0) then
         dzglu=z*glu 
         zlog = log((1d0+beta)/(1d0-beta))
         bgfcc = zlog*(z**2 + (1.d0-z)**2 +
     +        4d0*xksi*z*(1.d0-3.d0*z) - 8d0*xksi**2*z**2)
     +        + beta*(-1d0+8d0*z*(1d0-z)-z*(1d0-z)*4d0*xksi)
         bgfcc = 0.5d0*bgfcc*dzglu
	else
	   bgfcc=0d0
	endif
c	write(6,*) ' bgfcc: z ',bgfcc,z
	bgfcc = bgfcc*weight
c	write(6,*) ' bgfcc: x3 ,glu',bgfcc,x3,glu
c	write(6,*) ' bgfcc: x3,x,z ',x3,x,z
c	write(6,*) ' bgfcc: glu,dzglu ',glu,dzglu
      return
      end
      
	Subroutine nlogrid
      IMPLICIT None
      Double Precision X,xv,wmax
      DIMENSION X(20)
      Integer NDIM,NPOIN,NDIMEN,IMIX
      COMMON/DIVO/ NDIM,NPOIN
      COMMON /XVAL/ XV(20),NDIMEN
      REAL RQPMLO,RQPM,RQQB,RQQBC,RQQBB,RQCDC
      COMMON /OALPHAS/ WMAX,IMIX
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      double precision flux
      common/transform/flux
#include "rgqcdgrid.inc"
#include "rgraher.inc"
#include "rgpart.inc"
#include "rghs45.inc"
#include "rgrapgki.inc"
#include "rgpara.inc"
#include "rgpara1.inc"
#include "rgdisdif.inc"
#include "rgluco.inc"
#include "rglujets.inc"
#include "rgdiffr.inc"
#include "rgnlo.inc"
      double precision y,xxx,RSAFETY,xpo_max,xpo_min
      double precision RSUM,PQPMLO
	Double Precision dq2,testdx,dx
	Double Precision BGFNLO,BGFHQ
      Integer idiro,iproo,iheraco,nflavo,ndimeno,ndimo,ngo,npomo
      Integer nny,nnq,nnpom,ipom,IHERPYSO,iy,iq,i,iprint
	Double Precision ME,MP,W02,W12
	Real XPQ(-6:6),qtest
      DATA IPRINT/1/
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
	
	ME = P(1,5)
	MP = P(2,5)
	W02=(1.d0+MP)**2
      W12=W02-MP*MP
      YMAX=SSS+W12+DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMAX=YMAX/(2.D0*(SSS+ME*ME))
      YMIN=SSS+W12-DSQRT((SSS-W12)**2 - 4.D0*ME*ME*W12)
      YMIN=YMIN/(2.D0*(SSS+ME*ME))
      IF(YMI.GT.YMIN) YMIN=YMI
      IF(YMA.LT.YMAX) YMAX=YMA
 	write(6,*) ' nlogrid ymin,ymax,q2min,q2max',ymin,ymax,q2min,q2max 
     
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
	qtest = 10.
      IF(ID.EQ.1) THEN
         CALL rgstfu(2212,0.1,qtest,xpq)
         ELSE
         CALL RASTFU(KINT(2,2),0.1,qtest,xpq)
         ENDIF
      IF(IDISDIF.GE.1.AND.IDIRO.EQ.1) THEN
         IHERPYS = 0
         ID = 1
      ENDIF
      write(6,*) ' Now calcualting order NLO weights ',id
c      write(6,*) ' Please be patient, for diffraction it can take long'
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
         Q2MIN=ME*ME*MIN(YMIN*YMIN/(1.D0-YMIN),YMAX*YMAX/(1.D0-YMAX))
         IF(QMI.GT.Q2MIN) Q2MIN = QMI
         XPO_MIN=Q2MIN/YMAX/SSS
c 	   write(6,*) ' nlogrid xpom,xpmax,xpmin ',xpom,XPO_MAX,XPO_MIN
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
               Q2MIN=ME*ME*YY*YY/(1.D0-YY)
               IF(QMI.GT.Q2MIN) Q2MIN = QMI
               Q2MAX=YY*SSS - W12
               IF(QMA.LT.Q2MAX) Q2MAX = QMA
               Q2 = Q2MIN*(Q2MAX/Q2MIN)**X(2)      
		   flux=(4.d0*pi*alph**2)/(q2**2*y) 
     &         * ((y**2/2.d0) +1.0d0-y) 
               IF(ID.EQ.1) THEN
                  XPOM = 1.D0
               ELSE
                  XPOM = XPO_MIN*(XPO_MAX/XPO_MIN)**X(3)
               ENDIF
               xxx = q2/y/sss
c 	         write(6,*) ' nlogrid x,q2,y ',xxx,q2,y,sss
c 	      write(6,*) ' nlogrid xpom,xpmax,xpmin ',xpom,XPO_MAX,XPO_MIN
               IF(xxx.ge.1.) THEN
c                 write(6,*) 'iY,Iq',iy,iq
                  Q2=y*sss*0.99D0
                  goto 50
               endif
c Here get the total
               RQPMLO = SNGL(PQPMLO(xxx,Q2))
               RQQB=0.0
               RQQBC=0.0
               RQCDC=0.0
               if(RQPM.LE.0.0) THEN
                  RQPM = 1.e-10
                  GOTO 40
               ENDIF
c test collins calc.
		   dq2=dble(q2)
	         testdx=dq2/y/sss
	         dx =xxx	  
c	         write(6,*) ' bef bgf nlogrid x,q2 ',dx,dq2
 		
c	         print*,'dx',dx
c here you do the transformation to the cross section in nb/q**2
               RQQB=0d0
	         IF(dX.LT.1.and.dX.gt.0.) THEN
	            RQQB=SNGL(bgfnlo(dx,dq2))
c		      print*,'bgfnlo,x,q2', rqqb,dx,q2
		   endif

C........ gamma gluon fusion charm
               KPA = 4
cccc            IF(NFLAVO.GE.4.AND.IHF.EQ.0) THEN
               IF(NFLAVO.GE.4) THEN
                  RQQBC=0d0
	            IF(dX.LT.1.and.dX.gt.0.) THEN
	               RQQBC=SNGL(bgfhq(dx,dq2,4))
c		         print*,'bgfnlo,x,q2', rqqb,dx,q2
		       endif
                  IF(IHF.EQ.1.AND.IHFLA.EQ.5) RQQBC=0
                  RQQBB=0d0
                  IF(NFLAVO.GE.5) THEN
	               IF(dX.LT.1.and.dX.gt.0.) THEN
	                  RQQBB=SNGL(bgfhq(dx,dq2,5))
c		            print*,'bgfnlo,x,q2', rqqb,dx,q2
		          endif
                      IF(IHF.EQ.1.AND.IHFLA.EQ.4) RQQBB=0
                  ENDIF
               ENDIF
   40          CONTINUE
c            write(6,*) 'ccbar RQPM  RQQBC',RQPM,RQQBC
c           write(6,*) ' gamma gluon fusion c cbar'
               QY(IY+1) = DBLE(YY)
               QQ(IQ+1) = Q2
               QPOM(IPOM+1) = XPOM
c            write(6,*) ' y = ',yy,' Q2 = ',Q2
c            write(6,*) ' RQPM = ',RQPM,' RQQB = ',RQQB,'RQCDC = ',RQCDC,
c     +      'RQQBC = ',RQQBC
               RQPM = DBLE(RQQB+RQCDC+RQQBC+RQQBB+RQPMLO)
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
                  write(6,*) ' NLOGRID: RSUM > RQPM ',RSUM,RQPM
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
	
      double precision function pqpmlo(xi,q2)
      implicit none
c     calculates the F_2 LO
      double precision f2lo,x,q2,q,GEV2NB
      double precision  flux,xi
      common/transform/flux
      double precision beta,scale2
      common /spom/beta,scale2
	Double Precision dxbs(-6:6)
	Real xpq(-6:6)
	Integer I
      DATA GEV2NB/.3893857D+6/
      Integer ID
      double precision xpom,xold
      COMMON /PQPMID/xpom,ID
#include "rgpart.inc"	
#include "rgnlogrid.inc"
#include "rgnlodiff.inc"
	
	pqpmlo = 0d0
      x=xi
      IF(ID.EQ.1) THEN
        ELSE
        XOLD=X
        X=XOLD/XPOM
        IF(X.LT.0.9999.and.X.gt.0.) THEN
          ELSE
	    return
c          print*, 'PQPMLO x,xpom,beta unphys  ',xold,xpom,x
         ENDIf
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      q=dsqrt(q2)
	IF(q2.le.0) then
	  write(6,*) ' pqpmlo q2 < 0 ',q2
	  endif
	if(x.le.xpdf(1)) then
	   print*,' pqpmlo  x out of physical range',x,xpdf(1)
	endif
      if(x.lt.0.9999.and.x.gt.0.) then
         else
	   if(x.gt.0.9999) x=0.9999d0
c         print*,' Sigbgf 1.  x out of physical range',x3,rx3
      endif
	if(x.le.xpdf(1)) then
	   print*,' pqpmlo   x out of physical range',x,xpdf(1)
	endif
      IF(ID.EQ.1) THEN
         CALL rgstfu(2212,sngl(x),sngl(q2),xpq)
         ELSE
         CALL RASTFU(KINT(2,2),sngl(x),sngl(q2),xpq)
         ENDIF
      do i=-6,6
	dxbs(i)=dble(xpq(i))/diff_flux
	enddo
      f2lo=(4.d0/9.d0*(dxbs(-2)+dxbs(2))) +1./9.*(dxbs(-1)+
     +      dxbs(+1)+2*dxbs(3)) 
c     pqpmlo: lo contribution dsigma/dydq2 in [1/q**4]
      pqpmlo=f2lo*flux*GEV2NB*diff_flux
c	write(6,*) ' pqpmlo ',pqpmlo,flux
      return
      end

      double precision function bgfnlo(xi,q2)
      implicit none
      double precision x,q2,qq,GEV2NB
	double precision F2BGF
      double precision  flux,xi
      common/transform/flux
      double precision beta,scale2
      common /spom/beta,scale2
      DATA GEV2NB/.3893857D+6/
      Integer ID
      double precision xpom,xold
      COMMON /PQPMID/xpom,ID
#include "rgpart.inc"	
#include "rgnlogrid.inc"
#include "rgnlodiff.inc"
	Integer Pychge
	Integer I
	External Pychge
	External F2BGF
	
	bgfnlo = 0d0
      x=xi
      IF(ID.EQ.1) THEN
        ELSE
        XOLD=X
        X=XOLD/XPOM
        IF(X.LT.0.9999.and.X.gt.0.) THEN
          ELSE
c          print*, 'BGFNLO x,xpom,beta unphys  ',xold,xpom,x
	    return
         ENDIf
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      qq=dsqrt(q2)
	IF(q2.le.0) then
	  write(6,*) ' bgfnlo q2 < 0 ',q2
	  endif
	if(x.le.xpdf(1)) then
	   print*,' bgfnlo  x out of physical range',x,xpdf(1)
	endif
      if(x.lt.0.9999.and.x.gt.0.) then
         else
	   if(x.gt.0.9999) x=0.9999d0
c         print*,' bgfnlo 1.  x out of physical range',x3,rx3
      endif
	if(x.le.xpdf(1)) then
	   print*,' bgfnlo   x out of physical range',x,xpdf(1)
	endif
c       bgfnlo is the differential cross section  [1/q**4]
      qq=dsqrt(q2)
	bgfnlo = 0d0
	do i=-3,3
	if(i.ne.0) then
	   bgfnlo=bgfnlo+F2BGF(i,x,qq)*dble(pychge(i))**2
	endif
	enddo
	bgfnlo=bgfnlo/9.d0
	bgfnlo=flux*(bgfnlo)*GEV2NB*diff_flux
c	write(6,*) ' bgfnlo ',pqpmlo,flux
      return 
      end
          
      double precision function bgfhq(xi,q2,nfl)
      implicit none
      double precision x,q2,qq,GEV2NB
	double precision F2BGF
      double precision  flux,xi
      common/transform/flux
      double precision beta,scale2
      common /spom/beta,scale2
      DATA GEV2NB/.3893857D+6/
      Integer ID
      double precision xpom,xold
      COMMON /PQPMID/xpom,ID
#include "rgpart.inc"	
#include "rgnlogrid.inc"
#include "rgnlodiff.inc"
	Integer nfl,ifl
	Integer Pychge
	External Pychge
	External F2BGF
	
	bgfhq = 0d0
      x=xi
      IF(ID.EQ.1) THEN
        ELSE
        XOLD=X
        X=XOLD/XPOM
        IF(X.LT.0.9999.and.X.gt.0.) THEN
          ELSE
c          print*, 'BGFhq x,xpom,beta unphys  ',xold,xpom,x
           return
         ENDIf
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      qq=dsqrt(q2)
	IF(q2.le.0) then
	  write(6,*) ' bgfhq q2 < 0 ',q2
	  endif
	if(x.le.xpdf(1)) then
	   print*,' bgfhq  x out of physical range',x,xpdf(1)
	endif
      if(x.lt.0.9999.and.x.gt.0.) then
         else
	   if(x.gt.0.9999) x=0.9999d0
c         print*,' bgfhq 1.  x out of physical range',x3,rx3
      endif
	if(x.le.xpdf(1)) then
	   print*,' bgfhq   x out of physical range',x,xpdf(1)
	endif
c       bgfhq is the differential cross section  [1/q**4]
      qq=dsqrt(q2)
	bgfhq = 0d0
	ifl=nfl
	bgfhq=F2BGF(ifl,x,qq)*2d0*dble(pychge(ifl))**2/9.d0
	bgfhq=flux*bgfhq*GEV2NB*diff_flux
c	write(6,*) ' bgfhq ',ifl,pychge(ifl)
      return 
      end


      double precision function p_gqq(z)
      implicit none
      double precision z
      p_gqq=0.5d0*(1.d0-2.d0*z+2.d0*z**2)
      return 
      end

	Double Precision function ctnlo(q2,k2) 
	Implicit None
c     NLO subtraction cut function C(Q_1^2)
c     with Q_1^2=K2
c     1st: simple theta function C(Q_1^2)=\theta(Q^2-Q_1^2)
      Double Precision q2,k2
	Double Precision q2a,a
	Data a/0.5/
#include "rgnlo.inc"
	Icutsel = 0
	ctnlo = 0d0
	If(Icutsel.eq.0) Then
	   if(q2.gt.k2) then
	      ctnlo = 1d0
	   endif
	Elseif(Icutsel.eq.1) Then
	   q2a=q2/a
	   if(q2a.gt.k2) Then
	      ctnlo = 1d0 - a*k2/q2
         endif
	Else
	   write(6,*) ' Isel for NLO cut function not implemented ',
     &   Icutsel
	   Stop
	Endif      
	Return
	End

*****************************************************************************
********  Function NLOpdf (iptrn, x,q)   *************************************
******************************************************************************
      FUNCTION NLOpdf (IPRTN, X, Q)
C   returns x*f(x)
C   Given the parton distribution function in the array Upd in
C   COMMON / NLOPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C     July 29, 2001, Xiaomin
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Double precision NLOpdf

C
      PARAMETER (MXX = 200, MXQ = 30, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First, Qfirst
      Common 
     > / NLOPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX), UCC(MXPQX)
     > / NLOPar2 / Nx, Nt
     > / NLOXQrange / Qini, Qmax, Xmin, Xcr
C
      Logical NLO, INI_NLO
      Common /NLOpara/ NLO, INI_NLO,NSCHEME, NfMx, ifil

      Dimension Fq(M1), Df(M1)

      Data First, Qfirst/.true., .true./
      save First, Qfirst
      Integer Ncountx,Ncountq, Nmax
	Data Ncountx/0/
	Data Ncountq/0/
	Data Nmax/10/

C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         Ncountx = Ncountx + 1
         Print '(A, 2(1pE12.4))', 
     >    ' WARN_NLO: X << Xmin, extrapolation used; X, Xmin =', X, Xmin
         IF(Ncountx.gt.Nmax) Then
	     First = .false.
	     write(6,*) ' WARN_NLO: X << Xmin: Last Warning printed '
	     Endif
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif


C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif


      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini.and.Qfirst) then
            Ncountq = Ncountq + 1
            Print '(A, 2(1pE12.4))', 
     >     ' WARN_NLO: Q << Qini, extrapolation used; Q, Qini =', Q,Qini
            write(6,*) ' IPRTN, X, Q ',IPRTN, X, Q
c            Qfirst=.false.
         IF(Ncountq.gt.Nmax) Then
           Qfirst=.false.
	     write(6,*) ' WARN_NLO Q << Qini : Last Warning printed '
	     Endif
         endif
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >   ' WARN_NLO: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif


         Ip = Iprtn
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx 
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
c      write(6,*) 'Now interpolate in x for M1 Q"s'
	   Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue

C                                          Finish off by interpolating in Q
c      write(6,*) 'Finish off by interpolating in Q'
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)
      If ( Ftmp.lt.0.0d0) then
         NLOpdf=0.0d0
      else
         NLOpdf = Ftmp*x
      endif
      RETURN
      END

******************************************************************************
********  Function F2BGF(iptrn, x,q)   *************************************
******************************************************************************
      FUNCTION F2BGF (IPRTN, X, Q)
C
C   Given the parton distribution function in the array UCC in
C   COMMON / NLOPar1 / , this routine fetches deltabs(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C     July 29, 2001, Xiaomin
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Double precision F2BGF

C
      PARAMETER (MXX = 200, MXQ = 30, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First, Qfirst
      Common 
     > / NLOPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX), UCC(MXPQX)
     > / NLOPar2 / Nx, Nt
     > / NLOXQrange / Qini, Qmax, Xmin, Xcr
C
      Logical NLO, Ini_nlo
      Common /NLOpara/ NLO, INI_NLO,NSCHEME, NfMx, ifil

      Dimension Fq(M1), Df(M1)

      Data First /.true./
      Data Qfirst /.true./
      save First, Qfirst

c      if(NLO.and.(.not.Ini_nlo)) then
c         print *, 'nlopdf ', x, q, iprtn
c      endif
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
c         First = .false.
         Print '(A, 2(1pE12.4))', 
     >  ' WARN_F2BGF: X << Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif


C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif


      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini.and.Qfirst)  Then
            Qfirst=.false.
            Print '(A, 2(1pE12.4))', 
     > ' WARN_F2BGF: Q << Qini, extrapolation used; Q, Qini =', Q,Qini
         endif
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     > ' WARN_F2BGF: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif


c      If (Iprtn .GE. 3) Then
c         Ip = - Iprtn
c      Else
         Ip = Iprtn
c      EndIf
C                             Find the off-set in the linear array UCC
      JFL = Ip + NfMx 
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), UCC(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue

C                                          Finish off by interpolating in Q
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)
      If(Ftmp.lt.0.0d0) then
         F2BGF=0.0d0
      else
         F2BGF = Ftmp*x
      endif

      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
********************************************************************
C July 29, 2001, Xiaomin
C This subroutine does the following:
C     1. get msbar_pdf at Iparton, xx and QQ
C     2. transform  msbar_pdf to nlo_pdf at Iparton, x and Q 
C        where, x and Q are given by the grid we want.
C     3. Write the nlo_pdf in a array for later use.
C
C
********************************************************************
C  The organization of data in UPD(MXPQX) is the following.
C  from beginnig to end of the UPD array  !!!!This acctually correspond 
C  !!!!! to the acctual particle code of the event generator!!!!!!!!
C  1st block anti-top 
C  2nd block anti-botom 
C  3rd block anti-charm 
C  4th block anti-strane
C  5th block anti-d
C  6th block anti-u
C  7th block gluon
C  8th block u
C  9th block d
C 10th block s
C 11th block c
C 12th block b
C 13th block t
C..... Within each block the organization of the data is the following
c.....  Q=1  x=1E-5, x= 2E-5 .......  x=1
c.....  Q=2  x=1E-5, x= 2E-5 .......  x=1
c.....  Q=3  x=1E-5, x= 2E-5 .......  x=1
c.....  ................................
c.....Q=1E4  x=1E-5, x= 2E-5 .......  x=1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Nfmx=  # of flavor of quarks //// this has to be given by event generator.
c         gluon is only one.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine NLOset 

      Implicit Double Precision (A-H,O-Z)

      PARAMETER (MXX = 200, MXQ = 30, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)

      Common 
     > / NLOPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX), UCC(MXPQX)
     > / NLOPar2 / Nx, Nt
     > / NLOXQrange / Qini, Qmax, Xmin, Xcr

      Logical Ini_NLO, NLO
      Common /NLOpara/ NLO, INI_NLO,NSCHEME, NfMx, ifil
      Dimension UPDs(MXPQX), UCCs(MXPQX)

      Double precision NLOtrans
      Double precision NLOtmp 

      Data Nx, Nt,Al/ 90, 25, 0.25d0/
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c......find the grid value for Q, which is linear in ln(lnQ/Al)
c......ql(i) stores the value ln(q_i/al)
      delta=log(log(Qmax/Al))-log(log(Qini/AL))
      ql(0)= log(Qini/Al)+0.0001d0
      do 1 I=1, Nt
         y = dble(I) /dble(Nt) * delta
         ql(I) = exp(y) + ql(0)
         if(ql(I) .ge. log(Qmax/Al) ) ql(I) = log((qmax-1.0d0)/al)
 1    continue

c......find the grid value for x, follows Cteq convention
      do 5 I=0, Nx
         z = dble(I)/dble(Nx)
         xv(I) = XFRMZ1(z, Xcr, Xmin)
         if(xv(I) .ge. 1.0d0) xv(I) = 0.99999
         if(xv(I) .le. Xmin) xv(I) = Xmin
 5    continue


C             If data file not initialized, do so.
      If(ifil.eq.0) then


c........Fill the table UPD
         Do Iptr= -(Nfmx), (Nfmx)
            write(6,*) ' in NLOset Iptr loop, Nt,Nx  ', Iptr,Nt,Nx
            Do Iq = 0, Nt
               qq = exp(QL(Iq))*Al
c               write(6,*) 'in NLOset, qq', qq, al
               Do Ix = 0, Nx
                  xx = XV(Ix)
                  Itmp = (Iptr+Nfmx)*(1+Nt)*(1+Nx) + Iq*(1+Nx) + Ix + 1
                  UPD(Itmp) = NLOtrans(Iptr, xx, qq)/xx
                  Ucc(Itmp) = NLOtmp(Iptr,xx,qq)/xx-UPD(Itmp)
c			write(6,*) Ucc(Itmp),xx,qq,NLOtmp(Iptr,xx,qq)
c                  write(6,*) 'in NLOset, qq, xx ', qq, xx,Iptr
	            write(32) Itmp,UPD(Itmp),Ucc(Itmp)  
               enddo
            enddo
         enddo
      else
c........ Read the table UPD from the file
         Do Iptr= -(Nfmx), (Nfmx)
            Do  Iq = 0, Nt
               qq = exp(QL(Iq))*Al
c               write(6,*) 'in NLOset, qq', qq, al
               Do Ix = 0, Nx
                  xx = XV(Ix)
                  Itmp = (Iptr+Nfmx)*(1+Nt)*(1+Nx) + Iq*(1+Nx) + Ix + 1
	            read(32) Itmps,UPDs(Itmp),Uccs(Itmp)  
                  UPD(Itmp) = UPDs(Itmp)
                  Ucc(Itmp) = Uccs(Itmp)
               enddo
            enddo
         enddo
	

      endif
 100  format(1PE12.5)

      Return
      End


*****************************************************************************
      DOUBLE PRECISION FUNCTION XFRMZ1 (ZZ, XXCR, XXMIN)
*****************************************************************************
C                                 Invert the equation ZZ = ZFRMX (XX)
C                                 XXMIN and XXCR are the two parameters in
C                                 the transformation formula
      IMPLICIT NONE
      DOUBLE PRECISION ZZ, XXCR, XXMIN

      COMMON / XFRMZ_COM / Z, XCR, XMIN
      DOUBLE PRECISION Z, XCR, XMIN

      EXTERNAL ZFXL1, ZBRNT
      DOUBLE PRECISION ZFXL1, ZBRNT

      DOUBLE PRECISION RELATIVE_ERROR, ZLOW, ZHIGH

      DOUBLE PRECISION XLA, XLB, RESULT
      INTEGER IRET

C                        Mathematically, the full range of Z is allowed.
C                        Physically, Z < 1. is the only constraint.
C                        Note however,  Z < 0  corresponds to X < Xmin.
C                        We allow the range ZLOW < Z < 1.
      Z = ZZ
      XCR = XXCR
      XMIN = XXMIN
      RELATIVE_ERROR = 1E-3
      ZLOW = -10
      ZHIGH = 1.00002

      IF (Z .LE. ZHIGH .AND. Z .GT. ZLOW) THEN
          XLA = LOG (XMIN) * 1.5
          XLB = 0.00001
          RESULT = ZBRNT (ZFXL1, XLA, XLB, RELATIVE_ERROR, IRET)
          IF (IRET .NE. 0) THEN
             PRINT *, 'Error in solving equation in XFRMZ: ',
     >             'Return code from ZBRNT = ', IRET
          ENDIF
      Else
          PRINT *, 'Z out of range in XFRMZ: ',
     >           'Z=', Z, ', ZMIN=', ZLOW, ', ZHIGH=', ZHIGH
          RESULT = 0
      EndIf

      XFRMZ1 = EXP(RESULT)
      END

C==================================================================
C==================================================================
      DOUBLE PRECISION FUNCTION ZFXL1 (XL)
      DOUBLE PRECISION XL
      COMMON / XFRMZ_COM / Z, XCR, XMIN
      DOUBLE PRECISION X, Z, XCR, XMIN
      DOUBLE PRECISION ZFROMX
      ZFROMX(X) = (XCR-XMIN) * LOG (X/XMIN) + LOG (XCR/XMIN) * (X-XMIN)
      X = EXP(XL)
      ZFXL1 = ZFROMX(X) / ZFROMX(1D0)  - Z
      END

C==================================================================

      FUNCTION ZBRNT(FUNC, X1, X2, TOL, IRT)
C                                                   -=-=- zbrnt
 
C                          Return code  IRT = 1 : limits do not bracket a root;
C                                             2 : function call exceeds maximum
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)
      external func
 
      IRT = 0
      TOL = ABS(TOL)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        PRINT *, 'Root must be bracketed for ZBRNT. Set = 0'
        IRT = 1
        ZBRNT=0.
        RETURN
      ENDIF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRNT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      PRINT *, 'ZBRNT exceeding maximum iterations.'
      IRT = 2
      ZBRNT=B
      RETURN
      END


      
      
      
