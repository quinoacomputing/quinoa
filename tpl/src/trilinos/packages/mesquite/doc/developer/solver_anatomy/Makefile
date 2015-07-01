# Commands
LATEX = latex
DIA = dia
DVIPS = dvips
PS2PDF = ps2pdf13
EPSTOPDF = epstopdf
PDFLATEX = pdflatex
VIEWER = evince
SED = sed
TAR = tar
GZIP = gzip

# flags for PS2PDF
PDFFLAGS = -dEmbedAllFonts=true

# Document file name
DOCNAME = solver_anatomy

# Mesquite library for test-compiling example code
MESQUITE=../../clean/inst
MESQUITE_INC=$(MESQUITE)/include
MESQUITE_LIB=$(MESQUITE)/lib

# Some short-hand targets
.PHONEY: pdf ps view clean veryclean dist
pdf: $(DOCNAME).pdf
ps: $(DOCNAME).ps
view: $(DOCNAME).pdf
	$(VIEWER) $<
dist: $(DOCNAME).tgz

# clean
clean:
	rm -f *.aux \
              *.dvi \
              *.pdf \
              *.ps \
              *.lof \
              *.log \
              *.toc \
              *.tgz \
	      *.o \
	      figures.make \
	      $(FIGURES_PDF)

veryclean: clean
	rm -f $(FIGURES_EPS) figures.make

# tarball
$(DOCNAME).tgz: $(DOCNAME).tex $(FIGURES_EPS) figures.make Makefile
	test -d $(DOCNAME) || mkdir $(DOCNAME)
	test -d $(DOCNAME)/fig || mkdir $(DOCNAME)/fig
	cp $^ *.dia $(DOCNAME)/
	cp fig/*.dia $(DOCNAME)/fig/ 
	$(TAR) cf - $(DOCNAME) | $(GZIP) > $(DOCNAME).tgz 
	rm -rf $(DOCNAME)

# Suffix Rules
.SUFFIXES: .dia .eps .pdf .tex. .toc .dvi .ps

.dia.eps: 
	$(DIA) -e $@ -t eps $<

.eps.pdf:
	$(EPSTOPDF) $<

.tex.toc: 
	@echo "FIGURES=${FIGURES_PDF}"
	$(PDFLATEX) $<

.toc.pdf:
	$(PDFLATEX) $(<:.toc=.tex)

.toc.dvi: 
	$(LATEX) $(<:.toc=.tex)

.dvi.ps: 
	$(DVIPS) $< -o $@

# Scrape tex file for list of included figures
include figures.make
FIGURES_EPS = $(FIGURES:.x=.eps)
FIGURES_PDF = $(FIGURES_EPS:.eps=.pdf)
figures.make: $(DOCNAME).tex
	echo "FIGURES = \\" > $@ 
	$(SED) $(DOCNAME).tex \
           -e 's/\\psfig{\(.*,\)*figure=\([^,]\+\)\(,.*\)}/\2\\/p' \
           -e 's/\\includegraphics\(\[[^]]*\]\)\?{\(.*\)}/\2.x\\/p' \
           -e 'd' >> $@ || rm $@
	echo >> $@

# Output files depend onf igures
$(DOCNAME).pdf: ${FIGURES_PDF} 
$(DOCNAME).toc: ${FIGURES_PDF} 


# Test example code for syntax errors
CPPFLAGS = -I$(MESQUITE_INC)
check: example1.o example2.o

