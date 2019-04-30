C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c this function calculates the proton structure function at all x for
c q**2 < 10 GeV**2 using the parametrisation of Donnachie and Landshoff
c in DAMTP 93-23 and M/C-TH 93/11
C
C   RENAMED TO HSF2DL BY H.S.
c
c-----------------------------------------------------------------------
      FUNCTION HSF2DL(Q2,X)
      implicit double precision (a-h,o-z)
      c=0.219744
      b=0.278516
      x0=0.071348
      d=15.8769
      fm2=0.302408
      xi1=x*(1.+16./q2)
      xi2=x*(1.+1.7/q2)
      if (xi1.ge.1.) f1=0.
      if (xi1.lt.1.) f1=0.027/(xi1**0.0808)*(1.-xi1)**7*q2/(q2+6.25)
      if (xi2.ge.1.) f2=0.
      if (xi2.lt.1.) f2=2.0*c/(9.*xi2**0.0808)*(1.-xi2)**7*q2/(q2+1.)
      aa=1./(1./(2.038*c)-0.173)
      bb=0.489*b
      phi=q2/(q2+bb)
      phis=q2/(q2+aa)
      xi=x*(1.+0.28/q2)
      if (xi.ge.1.0) then
      HSF2DL=0.
      return
      else
      ht=d*x**2*(1.-xi)**2/(1.+q2/fm2)
      if(xi.lt.x0) then
      f3=10.*c/(9.*xi**0.0808)*phis
      s=f3+f2+f1+ht
      HSF2DL=S+B*XI**0.4525*PHI
      else
      ru=3.0*x0/(1.-x0)+0.4525
      rd=4.0*x0/(1.-x0)+0.4525
      buu=2.0/(x0**ru*((1.-x0)**3/0.4525-1./ru+3.0*x0/(1.+ru)
     &-3.0*x0**2/(2.+ru)+x0**3/(3.+ru))+1./ru-3.0/(1.+ru)
     &+3.0/(2.+ru)-1.0/(3.+ru))
      bu=buu*x0**(ru-0.4525)*(1.-x0)**3
      bdd=1.0/(x0**rd*((1.-x0)**4/0.4525-1./rd+4.0*x0/(1.+rd)
     &-6.0*x0**2/(2.+rd)+4.0*x0**3/(3.+rd)-x0**4/(4.+rd))+1/rd
     &-4.0/(1.+rd)+6.0/(2.+rd)-4.0/(3.+rd)+1.0/(4.+rd))
      bd=bdd*x0**(rd-0.4525)*(1.-x0)**4
      u=buu*xi**ru*(1.-xi)**3*phi
      d=bdd*xi**rd*(1.-xi)**4*phi
      rs=9.0*x0/(1.-x0)+0.4525
      rss=7.*x0/(1.-x0)-0.0808
      bss=(b-4.*bu/9.-bd/9.)*x0**(0.4525-rs)/(1.-x0)**9
      css=c/(x0**(0.0808+rss)*(1.-x0)**7)
      f3=(bss*xi**rs*(1.-xi)**2*phi+(10./9.)*css*xi**rss*phis)
     &*(1.-xi)**7
      s=f3+f2+f1+ht
      HSF2DL=4.*U/9.+D/9.+S
      return
      endif
      endif
      return
      end
