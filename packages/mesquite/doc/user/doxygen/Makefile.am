USER_GUIDE_PDF = @USER_GUIDE_PDF@
USER_GUIDE_PS = @USER_GUIDE_PS@

if WITH_DOXYGEN_PDF
  pdf_DATA = MesquiteAPI.pdf
else
  pdf_DATA =
endif

if WITH_DOXYGEN_PS
  ps_DATA = MesquiteAPI.ps
else
  ps_DATA =
endif

htmldir = @htmldir@/html
if WITH_DOXYGEN_HTML
  html_DATA = html/*
else
  html_DATA =
endif

.PHONY: html html/*
html: html/index.html
html/*: html/index.html

html/index.html latex/refman.tex: Mesquite-user.dox
	$(DOXYGEN) $<

latex/refman.dvi: latex/refman.tex latex/doxygen.sty
	cd latex && $(LATEX) refman.tex	
	cd latex && makeindex refman.tex
	cd latex && $(LATEX) refman.tex	
	cd latex && $(LATEX) refman.tex	
	cd latex && $(LATEX) refman.tex	
	cd latex && $(LATEX) refman.tex	

latex/refman.ps: latex/refman.dvi
	cd latex && $(DVIPS) refman.dvi

MesquiteAPI.ps: latex/refman.ps
	cp $< $@

MesquiteAPI.pdf: MesquiteAPI.ps
	ps2pdf $< $@


install-pdf-local: install-pdfDATA

install-ps-local: install-psDATA

install-dvi-local: install-dviDATA

install-html-local: install-htmlDATA html
