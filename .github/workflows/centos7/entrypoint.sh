#!/bin/sh -l
set -x
uname -a 
cat /etc/issue
yum -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
yum -y install yum-plugin-copr
yum -y copr enable averbyts/fastjet
yum -y install  git  pythia6 
yum -y install  HepMC HepMC-devel
yum -y install  HepMC3 HepMC3-devel HepMC3-search HepMC3-search-devel
yum -y install  Rivet Rivet-devel 
yum -y install  lhapdf lhapdf-devel
yum -y install  gcc gcc-c++  gcc-gfortran make
yum -y install  autoconf automake libtool
yum -y install  zlib zlib-devel
yum -y install  texlive-latex-bin-bin texlive-metafont-bin ghostscript texlive-fontspec texinfo 


cd rapgap-3.303
#No idea why manual buid fails on CI
sed -i 's/manual//g' Makefile.am

autoreconf -fisv
sed -i 's/[^[:print:]\r\t]//g' src/*/*f
sed -i 's/[^[:print:]\r\t]//g' src/*/*F
sed -i 's/[^[:print:]\r\t]//g' bases51/*f   
sed -i 's/[^[:print:]\r\t]//g' misc/*F   

./configure --with-lhapdf=/usr --with-hepmc=/usr  --with-pythia6=/usr --with-rivet=/usr --with-hepmc3=/usr
make

src/rapgap_hepmc  < data/steer-ep-no-qedrad-dis-rivethepmc23
head -n 30 hm3output.hepmc3
out=$?
echo ::set-output name=out::$out
