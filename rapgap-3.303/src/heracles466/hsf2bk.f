
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      SUBROUTINE HSF2BK(MODE,X,Q2,F2P,F2D)
C     SUBROUTINE F2PD(MODE,X,Q2,F2P,F2D)
c
c --- A code to calculate the nucleon structure function F2 in the
c --- low q2 and low x region within the GVMD inspired model described in
c --- B.Badelek and J.Kwiecinski, B295 (1992) 263.
c --- The structure functions F2 corresponding to the QCD improved
c --- parton model are obtained from the interpolation formula based on
c --- Tchebyshev polynomials. The vector meson contributions are calculated
c --- directly. Parton model contributions are based on the MRS
c --- parton distributions, A.D.Martin, W.J.Stirling and R.G.Roberts,
c --- PR D47(1993) 145, Phys.Lett B306(1993) 145.
c
c >>>>>>  Model is valid for
c                             nu > 10 GeV,
c                             10**(-5) < x < 0.1,
c                             Q2< 1000 GeV2                  <<<<<<
c
c --- The following data files exist:
c --- mode=1. f2mod1.dat, corresponds to the unshadowed D- parametrisation;
c --- mode=2. f2mod2.dat, to the shadowed D- parametrisation and R=2GeV**(-1);
c --- mode=3. f2mod3.dat, to the shadowed D- parametrisation and R=5GeV**(-1);
c --- mode=4. f2mod4.dat, to the unshadowed D0 set of partons;
c --- mode=5. f2mod5.dat, to the D-';
c --- modes 1-4 have calculations done in the LO and mode=5 in the NLO of QCD.
c --- OBS 1:  mode=1,2,3,4 have a narrower validity range: 10**(-4) < x < 0.1
c --- OBS 2:  x validity range applies in fact to the {\bar x} variable,
c             cf.the publication for the definition; validity in the Bjorken x
c             extends to even lower x values (e.g. to x=0 for photoproduction).
c --- User has to attach a data file containing the interpolation
c --- coefficients for a chosen parametrisation.
c
c --- Subroutine parameters:
c --- mode: an integer defining the input file (e.g.mode=2 for file f2mod2.dat),
c --- q2,x: the usual kinematic variables,
c --- f2p=F2 for the proton and f2d=F2 for the deuteron (returned values).
c --- Example of the usage of the subroutine is given above.
c
c --- History:
c     Oct. 25th, 1993; introduced mode=5 which corresponds to the MRS D-'
c                      parametrisation; NLO evolution (LO before);
c                      extension of model validity down to x=1.E-5
c                      (old limit: x=1.E-4)
c
      common//h(2,20,20)
      dimension alama(5),ftwo(2)
      data nmax, pmass / 20, 0.938272/
      data w02, ymax, q20, q2f / 1.2, 9.2103, 1., 1000./
      data alama/3*0.2304, 0.1732, 0.2304/
c
      qbar=q2+w02
      cbar=q2/qbar
      viu=q2/x
      viug=viu/(2.*pmass)
      viubar=viu+w02
      xbar=qbar/viubar
c
      alam=alama(mode)*alama(mode)
      ql0=alog(q20/alam)
      qlf=alog(q2f/alam)
      qlp=ql0+qlf
      qlm=qlf-ql0
      aq=acos((2.*alog(qbar/alam)-qlp)/qlm)
      if (mode.gt.4) ymax=11.5129
      ax=acos((2.*alog(1./xbar)-ymax)/ymax)
c
c
      do 500 j=1,2
      f2=0.
      do 400 nx=1,nmax
      anx=float(nx)-1.
      tx=cos(anx*ax)
      do 300 nq=1,nmax
      anq=float(nq)-1.
      tq=cos(anq*aq)
      f2=f2+tx*tq*h(j,nq,nx)
300   continue
400   continue
      ftwo(j)=4./float(nmax)/float(nmax)*f2
500   continue
c
c --- vector meson contribution
c
      call sigvmes(viug,sigro,sigfi)
      call vmesnuc(q2,sigro,sigfi,fv)
c
      f2p=cbar*ftwo(1)+fv
      f2d=cbar*ftwo(2)+fv
c
      return
      end
