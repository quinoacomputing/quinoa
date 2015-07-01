if ENABLE_DOXYGEN
  DOXY_SUBDIR = doxygen
endif
SUBDIRS = $(DOXY_SUBDIR)

UG_SRCS = users-guide.tex \
	ack.tex \
	api.tex \
	caveats.tex \
	custom.tex \
	data.tex \
	extend.tex \
	examples.tex \
	geom.tex \
	install.tex \
	intro.tex \
	mesh.tex \
	more_tmp.tex \
	opt_strategies.tex \
	plan.tex \
        parallel.tex \
	debugging.tex \
	support.tex \
	team.tex \
	tmp.tex \
        trouble.tex \
	wrappers.tex \
        MesquiteUI.eps \
        hole_in_square.ps \
        msq-paradigm.eps \
        square_rand.eps \
        parallel_mesh.eps \
	generalized_optimization.eps \
	patches-horiz.eps \
	iterplot.eps \
	mesh_gpt.eps \
	mesh_vtk.eps \
	mesh.eps \
	meshqual.eps \
	meshqual3d.eps \
        mesquite.bib

EXTRA_DIST = $(UG_SRCS) \
	     README.USERS-GUIDE

if WITH_USER_GUIDE_PDF
  pdf_DATA = users-guide.pdf
endif
if WITH_USER_GUIDE_PS
  ps_DATA = users-guide.ps
endif

SUFFIXES = .dvi .pdf .ps
.ps.pdf:
	$(PS2PDF) $< $@

.dvi.ps:
	$(DVIPS) $< -o $@

users-guide.dvi: users-guide.tex $(UG_SRCS) 
	$(LATEX) $<
	$(BIBTEX) users-guide.aux
	$(LATEX) $<
	$(LATEX) $<

pdf: users-guide.pdf

ps: users-guide.ps


install-pdf-local: install-pdfDATA

install-ps-local: install-psDATA

install-dvi-local: install-dviDATA

MOSTLYCLEANFILES = *.aux *.bbl *.blg *.lof *.log *.lot *.toc 
CLEANFILES = users-guide.dvi users-guide.pdf users-guide.ps


