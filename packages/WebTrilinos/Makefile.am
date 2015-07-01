# @HEADER
# ************************************************************************
#
#                WebTrilinos: A Web Interface to Trilinos
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
#
# ************************************************************************
# @HEADER

## #######################################################################
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Aclocal command (rarely used - don't worry about it)
## #######################################################################
ACLOCAL_AMFLAGS = -I config

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = \
config/depcomp            config/install-sh \
config/missing            config/mkinstalldirs \
config/generate-makeoptions.pl \
config/replace-install-prefix.pl \
config/config.guess       config/config.sub \
bootstrap \
index.html                config/token-replace.pl \
ChangeLog-WebTrilinos     README-WebTrilinos \
common/index.html         common/basic.html \
common/setup.html         common/testing.html \
common/favicon.ico        common/footer.html \
common/util.php           common/config.php \
common/load_images.js     common/header.html \
common/header.php         common/styles.css \
common/WebTrilinos_config.h.in \
data/index.html \
data/HB/index.html \
data/HB/fidap005.rua \
data/XML/index.html \
data/XML/sample.xml \
data/MM/index.html \
tmp/index.html \
img/beginner.png             img/rd100.png \
img/trilinos_normal.png      img/recirc1.png \
img/trilinos_over.png        img/eth_logo.gif \
img/recirc2.png              img/trilinos_title_normal.png \
img/expert.png               img/sandia_normal.png \
img/trilinos_title_over.png  img/intermediate.png \
img/sandia_over.png          img/webtrilinos_logo.jpg \
img/laplace2d_stencil.png    img/spacer.png \
img/matrix_portal_logo.png   img/star2d_stencil.png \
MatrixPortal/calculator.html MatrixPortal/index.html \
MatrixPortal/step1.html      MatrixPortal/step1C.html \
MatrixPortal/step2.html      MatrixPortal/step2.html \
MatrixPortal/step3.html      MatrixPortal/step4.html \
MatrixPortal/MATLAB_to_XML.m \
MatrixPortal/step_process.py \
MatrixPortal/step_test.txt \
MatrixPortal/help/footer.html \
MatrixPortal/help/header.html \
MatrixPortal/help/step1_collections.html \
MatrixPortal/help/step1_generator.html \
MatrixPortal/help/step1_upload.html \
MatrixPortal/help/step2_multilevel.html \
MatrixPortal/help/step2_incomplete_factorization.html \
MatrixPortal/help/step2_performance.html \
MatrixPortal/help/step2_general.html \
MatrixPortal/help/step2_relaxation.html \
MatrixPortal/help/step2_chebyshev.html \
MatrixPortal/help/step4_plot.html \
MatrixPortal/help/step_workflow.html \
MatrixPortal/help/template.html

## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs 
#
#  Again, I hope that AUX_DIR(../../config) eliminates these
#  config/install-sh config/missing config/mkinstalldirs 

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	$(AUX_DIST) 

all-local:
	@echo "        +----------------------------------------------------+"
	@echo "        | Trilinos package WebTrilinos built successfully    |"
	@echo "        |                                                    |"
	@echo "        | *Note that WebTrilinos must be manually installed* |"
	@echo "        | *to be operational through a web browser*          |"
	@echo "        |                                                    |"
	@echo "        | Details are reported on the web site               |"
	@echo "        | http://trilinos.sandia.gov/packages/               |"
	@echo "        |        webtrilinos/installation.html               |"
	@echo "        |                                                    |"
	@echo "        | Please send questions/comments to:                 |"
	@echo "        | - WebTrilinos-users@software.sandia.gov            |"
	@echo "        | - WebTrilinos-developers@software.sandia.gov       |"
	@echo "        |                                                    |"
	@echo "        | Have fun!                                          |"
if HAVE_MPI
	@echo "        |                                                    |"
	@echo "        | NOTE: YOU HAVE CONFIGURED TRILINOS WITH SUPPORT    |"
	@echo "        |       FOR MPI. WebTrilinos CURRENTLY DOES NOT      |"
	@echo "        |       SUPPORT MPI. YOU WILL BE ABLE TO COMPILE AND |"
	@echo "        |       RUN ALL THE EXAMPLES, BUT INSTALLING         |"
	@echo "        |       WebTrilinos ON A WEB SERVER WILL REQUIRE TO  |"
	@echo "        |       HACK FEW MPI COMMANDS (LIKE lamboot/lamhalt) |"
	@echo "        |       IN THE SCRIPTS. SEE DETAILS ON THE WEB PAGES |"
	@echo "        |       OR RE-CONFIGURE *WITHOUT* MPI SUPPORT.       |"
endif
	@echo "        +----------------------------------------------------+"
	@echo 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though BUILD_TESTS and BUILD_EXAMPLES will never be
#defined, so that the tests and examples are included in the distribution
#tarball.

if SUB_TEST
TEST_SUBDIR=test c++
endif

if SUB_EXAMPLE
EXAMPLE_SUBDIR=python
endif

if BUILD_TESTS
tests:
	@echo ""
	@echo "Now building WebTrilinos tests."
	@echo ""
	cd $(top_builddir)/test && $(MAKE) && cd ../c++ && $(MAKE)
	@echo ""
	@echo "Finished building WebTrilinos tests."
	@echo ""
else
tests:
	@echo "WebTrilinos tests were disabled at configure time"
endif

if BUILD_EXAMPLES
examples:
	@echo ""
	@echo "Now building WebTrilinos examples."
	@echo ""
	cd $(top_builddir)/python && $(MAKE) && cd .. && $(MAKE) tests
	@echo ""
	@echo "Finished building WebTrilinos examples."
	@echo ""

else
examples:
	@echo "WebTrilinos examples were disabled at configure time"
endif

install-examples:
	@echo "No WebTrilinos examples should be installed"

clean-tests:
	cd $(top_builddir)/test && $(MAKE) clean && cd ../c++ && $(MAKE) clean

clean-examples:
	cd $(top_builddir)/python && $(MAKE) clean && cd .. && $(MAKE) clean-tests

everything:
	$(MAKE) && $(MAKE) examples && $(MAKE) tests

clean-everything:
	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

install-everything:
	$(MAKE) install && $(MAKE) install-examples

SUBDIRS = $(EXAMPLE_SUBDIR) $(TEST_SUBDIR)

## #######################################################################
## Other make targets
## #######################################################################

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..

if HAVE_MPI
WEBTRILINOS_CHECK_COMM=mpi
else
WEBTRILINOS_CHECK_COMM=serial
endif

run-tests:
	@echo
	@echo "================================"
	@echo "Running the WebTrilinos tests..."
	@echo "================================"
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(WEBTRILINOS_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=WebTrilinosTests \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=WebTrilinos

run-pytests:
	@echo
	@echo "======================================="
	@echo "Running the Python/WebTrilinos tests..."
	@echo "======================================="
	@echo
	@echo "For more details, please refer to Trilinos/doc/check_installation.txt"
	@echo
	@$(PERL_EXE) $(TRILINOS_HOME_DIR)/commonTools/test/utilities/runtests \
	--trilinos-dir=$(TRILINOS_HOME_DIR) \
	--comm=$(WEBTRILINOS_CHECK_COMM) \
	--mpi-go="$(TRILINOS_MPI_GO)" \
	--build-dir=$(TRILINOS_BUILD_DIR) \
	--category=WebTrilinosPyTests \
	--output-dir=@abs_top_builddir@/test/runtests-results \
	--verbosity=1 \
	--packages=WebTrilinos
