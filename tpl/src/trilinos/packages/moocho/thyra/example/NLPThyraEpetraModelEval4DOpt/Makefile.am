# @HEADER
# ************************************************************************
# 
# Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
#                  Copyright (2003) Sandia Corporation
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
# Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
# 
# ************************************************************************
# @HEADER

all-local: test_lumped_param_model.py moochoThyraFdNand.xml

#
# Build the external example for users to see how to do this.
#

if USING_EXPORT_MAKEFILES
EXTERNAL_MAKEFILE_EXAMPLE = external
endif

SUBDIRS = $(EXTERNAL_MAKEFILE_EXAMPLE)

#
# Build the actual test program in the build directory so that
# we can test this.
#

if HAVE_MPI
NLPTHYRAEPETRAMULTIPOINTMODELEVAL4DOPT = NLPThyraEpetraMultiPointModelEval4DOpt
endif

noinst_PROGRAMS = NLPThyraEpetraModelEval4DOpt $(NLPTHYRAEPETRAMULTIPOINTMODELEVAL4DOPT)

NLPThyraEpetraModelEval4DOpt_SOURCES = \
  $(srcdir)/NLPThyraEpetraModelEval4DOptMain.cpp\
  $(top_srcdir)/../epetraext/example/model_evaluator/4dopt/EpetraModelEval4DOpt.cpp
NLPThyraEpetraModelEval4DOpt_DEPENDENCIES = $(common_dependencies)
NLPThyraEpetraModelEval4DOpt_LDADD = $(common_ldadd)

NLPThyraEpetraMultiPointModelEval4DOpt_SOURCES = \
  $(srcdir)/NLPThyraEpetraMultiPointModelEval4DOptMain.cpp\
  $(top_srcdir)/../epetraext/example/model_evaluator/4dopt/EpetraMultiPointModelEval4DOpt.cpp
NLPThyraEpetraMultiPointModelEval4DOpt_DEPENDENCIES = $(common_dependencies)
NLPThyraEpetraMultiPointModelEval4DOpt_LDADD = $(common_ldadd)

test_lumped_param_model.py: $(builddir)/_test_lumped_param_model.py.stub
	cp $(builddir)/_test_lumped_param_model.py.stub test_lumped_param_model.py
	perl -i -pe 's/MOOCHO_NLPThyraEpetraModelEval4DOpt/NLPThyraEpetraModelEval4DOpt/g' test_lumped_param_model.py
	chmod a+x test_lumped_param_model.py

moochoThyraFdNand.xml: $(srcdir)/_moochoThyraFdNand.xml
	cp $(srcdir)/_moochoThyraFdNand.xml moochoThyraFdNand.xml

include $(top_srcdir)/MakefileExeStub.mak

AM_CPPFLAGS += -I$(top_srcdir)/../epetraext/example/model_evaluator/4dopt


# Extra files to add to the tarball
EXTRA_DIST = _test_lumped_param_model.py.stub.in _moochoThyraFdNand.xml
