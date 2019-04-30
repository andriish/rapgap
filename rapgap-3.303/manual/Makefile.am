#if BUILD_DOCS
pkgdocdir = ${datadir}
#dist_
pkgdoc_DATA = rapgap32.pdf
EXTRA_DIST = rapgap32.tex  mcite.sty dina4.sty dismc.epsi difmc.epsi
CLEANFILES =rapgap32.pdf
DISTCLEANFILES = rapgap32.log rapgap32.idx rapgap32.aux rapgap32.dvi \
                 rapgap32.ilg rapgap32.blg rapgap32.bbl rapgap32.ind copy-doc

all-local: copy-doc rapgap32.pdf

copy-doc:
	find ${srcdir} -type f -exec ln -s {} 2> /dev/null \;
	if test -h rapgap32.ind; then rm -f rapgap32.ind; cp ${srcdir}/rapgap32.ind . && chmod a+w rapgap32.ind; fi
	touch copy-doc

rapgap32.pdf: copy-doc
	latex rapgap32  && latex rapgap32 && makeindex rapgap32 && latex rapgap32 && dvipdf rapgap32


#endif
