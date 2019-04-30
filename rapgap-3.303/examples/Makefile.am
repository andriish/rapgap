## Copyright (C) 2006 Free Software Foundation, Inc.
## This Makefile.am is free software; the Free Software Foundation
## gives unlimited permission to copy, distribute and modify it.

AM_FFLAGS = -I$(top_srcdir)/include -I$(PYTHIA_SOURCE) $(MYFFLAGS)
AM_CPPFLAGS = -I$(srcdir)/../include -I$(HEPMC)/include
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = rapgap_hepmc

rapgap_hepmc_SOURCES =  rgmain.F   \
          			$(top_srcdir)/src/rapgap/upinit.F   $(top_srcdir)/src/rapgap/upevnt.F
	
LDADD = -L../src -lrapgap33 \
        -L../bases51 -lbases \
        $(HEPMCLIB) $(ARIADNELIB) $(PYTHIALIB) $(PDFLIB) ../misc/libmycern.la -lstdc++ -lm
        
all:
#	cp rapgap_hepmc ../bin/.        
