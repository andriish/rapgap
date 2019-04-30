C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      subroutine sigvmes(viu,sigro,sigfi)
c
c-----calculates sigro, sigfi at given energy viu using the energy dependent
c-----piN and KN total cross-sections
c
      data      amk,     ampi,pmbnag,sigpi0, sigk0,   czad,sigpi1,sigpi2
     *    /0.493646, 0.139567,  2.56, 13.52, 14.17, 0.0102, 22.77, 2.35/
      data sigk1, alfpi1, alfpi2, alfk1
     *    /14.75,  0.369,   0.37, 0.515/
c
      gampi=viu/ampi
      gamk=viu/amk
      gapa=alog(gampi)
      gaka=alog(gamk)
      spip=sigpi1*viu**(-alfpi1)
      spip=spip+sigpi0*(1.+czad*gapa*gapa)
      spim=sigpi2*viu**(-alfpi2)
      skp=sigk1*viu**(-alfk1)
      skp=skp+sigk0*(1.+czad*gaka*gaka)
      sigro=pmbnag*spip
      sigfi=pmbnag*(2.*skp-spip-spim)
      return
      end
