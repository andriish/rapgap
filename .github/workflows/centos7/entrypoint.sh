#!/bin/sh -l
set -x
uname -a 
cat /etc/issue
yum -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
yum -y install yum-plugin-copr
yum -y copr enable averbyts/fastjet
yum -y install  git  pythia6 
yum -y install  HepMC HepMC-devel
yum -y install  HepMC3 HepMC3-devel
yum -y install  Rivet Rivet-devel 
yum -y install  lhapdf lhapdf-devel
yum -y install  gcc gcc-c++ 


autoreconf -fisv
sed -i 's/[^[:print:]\r\t]//g' src/*/*f
sed -i 's/[^[:print:]\r\t]//g' src/*/*F
sed -i 's/[^[:print:]\r\t]//g' bases51/*f   
sed -i 's/[^[:print:]\r\t]//g' misc/*F   

./configure --with-lhapdf=/usr --with-hepmc=/usr  --with-pythia6=/usr --with-rivet=/usr --with-hepmc3=/usr
make

out=$?
echo ::set-output name=out::$out
