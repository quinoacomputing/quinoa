#@HEADER
# ************************************************************************
# 
#               Didasko : a Tutorial for Trilinos
#                 Copyright (2005) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.didasko

EXAMPLE = ex1 ex2 ex3
EXEEXT = .exe

noinst_PROGRAMS = $(EXAMPLE)

ex1_SOURCES = $(srcdir)/ex1.cpp
ex1_CXXFLAGS = $(EXPORT_INCLUDES) 
ex1_LDADD = $(EXPORT_LIBS)

ex2_SOURCES = $(srcdir)/ex2.cpp
ex2_CXXFLAGS = $(EXPORT_INCLUDES) 
ex2_LDADD = $(EXPORT_LIBS)

ex3_SOURCES = $(srcdir)/ex3.cpp
ex3_CXXFLAGS = $(EXPORT_INCLUDES) 
ex3_LDADD = $(EXPORT_LIBS)

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(DIDASKO_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(DIDASKO_LIBS))
else
EXPORT_INCLUDES = $(DIDASKO_INCLUDES)
EXPORT_LIBS = $(DIDASKO_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

