

#
# A) Define the package
#

TRIBITS_PACKAGE(Didasko)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Latex
  HAVE_DIDASKO_LATEX
  "Enables compilation of LaTeX source files.[Note, this is currently unsupported]"
  OFF )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  examples/teuchos/xml-data
  src/.*.eps
  src/.*.tex
  src/.*.ps
  src/DOEbwlogo.pdf
  src/Makefile
  src/SANDreport.cls
  src/Trilinos60Tutorial.pdf
  src/Trilinos70Tutorial.pdf
  src/Trilinos80Tutorial.pdf
  src/TrilinosTutorial_ReviewAndApproval.doc
  src/chapterbox.pdf
  src/colabarticle.cls
  src/snllineblk.pdf
  src/tutorial_biblio.bib
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
