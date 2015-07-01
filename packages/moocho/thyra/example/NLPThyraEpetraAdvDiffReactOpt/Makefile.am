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

all-local: square.1.000 _run-varied-tests.pl aztecooParams.xml

noinst_PROGRAMS = \
  NLPThyraEpetraAdvDiffReactOpt \
  MultiPeriodNLPThyraEpetraAdvDiffReactOpt

NLPThyraEpetraAdvDiffReactOpt_SOURCES = \
  $(common_sources) \
  $(srcdir)/NLPThyraEpetraAdvDiffReactOptMain.cpp

NLPThyraEpetraAdvDiffReactOpt_DEPENDENCIES = $(common_dependencies)
NLPThyraEpetraAdvDiffReactOpt_LDADD = $(common_ldadd)

MultiPeriodNLPThyraEpetraAdvDiffReactOpt_SOURCES = \
  $(common_sources) \
  $(srcdir)/MultiPeriodNLPThyraEpetraAdvDiffReactOptMain.cpp

MultiPeriodNLPThyraEpetraAdvDiffReactOpt_DEPENDENCIES = $(common_dependencies) square.1.000
MultiPeriodNLPThyraEpetraAdvDiffReactOpt_LDADD = $(common_ldadd)

square.1.000: $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/data/square/square.1.000
	cp $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/data/square/square.1.000 square.1.000
	cp $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/data/square/square.2.000 square.2.000
	cp $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/data/square/square.9x9.2.000 square.9x9.2.000
	cp $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/data/square/square.9x9.2.001 square.9x9.2.001

_run-varied-tests.pl : run-varied-tests.pl
	cp $(srcdir)/run-varied-tests.pl _run-varied-tests.pl

aztecooParams.xml: _aztecooParams.xml
	cp $(srcdir)/_aztecooParams.xml aztecooParams.xml

include $(top_srcdir)/MakefileExeStub.mak

CPPFLAGS += -I$(top_srcdir)/../thyra/example/operator_vector \
  -I$(top_srcdir)/../epetraext/example/model_evaluator/GLpApp

common_sources = \
  $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/GLpApp_AdvDiffReactOptModel.cpp\
  $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/GLpApp_AdvDiffReactOptModelCreator.cpp\
  $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/GLpApp_GLpYUEpetraDataPool.cpp\
  $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/GenSQP_YUEpetraVector.cpp\
  $(top_srcdir)/../epetraext/example/model_evaluator/GLpApp/rect2DMeshGenerator.cpp

EXTRA_DIST = run-varied-tests.pl _aztecooParams.xml
