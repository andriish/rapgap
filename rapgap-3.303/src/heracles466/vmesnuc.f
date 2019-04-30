C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      subroutine vmesnuc(q2,sigro,sigfi,fvmes2)
c
c-----calculates the vector meson contribution to the nucleon (!) f2.
c
      dimension vmass(3),vcoupl(3),sigin(3)
      data pi/3.1415926/, vmass/0.7683, 0.78195, 1.019412/
      data vcoupl/1.98, 21.07, 13.83/
c
      sigin(1)=sigro
      sigin(2)=sigro
      sigin(3)=sigfi
      sumv=0.
      do 100 l=1,3
      sv=vmass(l)*vmass(l)
      denv=pi*vcoupl(l)*(q2+sv)*(q2+sv)
      sumv=sumv+sv*sv*sigin(l)/denv
100   continue
      fvmes2=q2/4./pi*sumv
      return
      end
