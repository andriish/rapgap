
## Copyright (C) 2006 Free Software Foundation, Inc.
## This Makefile.am is free software; the Free Software Foundation
## gives unlimited permission to copy, distribute and modify it.

ACLOCAL_AMFLAGS = -I m4  

if HZTOOL_OPT
   MAYBE_OPT = examples
else
   MAYBE_OPT = 
endif
if ARI_OPT
   ARIADNE_OPT = ariadne
else
   ARIADNE_OPT = ariadne-dummy
endif
SUBDIRS = include misc bases51  $(ARIADNE_OPT) src $(MAYBE_OPT) data manual
#SUBDIRS = include misc bases51 src $(MAYBE_OPT)  data manual
# SUBDIRS = manual
#dist_doc_DATA = README
EXTRA_DIST = README ChangeLog
