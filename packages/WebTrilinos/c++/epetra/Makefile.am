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
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
#
# ************************************************************************
# @HEADER

include $(top_builddir)/Makefile.export

EXAMPLE = ex_Comm \
	  ex_Map \
	  ex_Vector \
	  ex_MultiVector \
	  ex_CrsMatrix \
	  ex_FECrsMatrix \
	  ex_VbrMatrix
EXEEXT = .exe

noinst_PROGRAMS = $(EXAMPLE)

ex_Comm_SOURCES = $(srcdir)/ex_Comm.cpp
ex_Comm_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_Comm_LDADD = $(WEBTRILINOS_LIBS)

ex_Map_SOURCES = $(srcdir)/ex_Map.cpp
ex_Map_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_Map_LDADD = $(WEBTRILINOS_LIBS)

ex_Vector_SOURCES = $(srcdir)/ex_Vector.cpp
ex_Vector_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_Vector_LDADD = $(WEBTRILINOS_LIBS)

ex_MultiVector_SOURCES = $(srcdir)/ex_MultiVector.cpp
ex_MultiVector_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_MultiVector_LDADD = $(WEBTRILINOS_LIBS)

ex_CrsMatrix_SOURCES = $(srcdir)/ex_CrsMatrix.cpp
ex_CrsMatrix_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_CrsMatrix_LDADD = $(WEBTRILINOS_LIBS)

ex_FECrsMatrix_SOURCES = $(srcdir)/ex_FECrsMatrix.cpp
ex_FECrsMatrix_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_FECrsMatrix_LDADD = $(WEBTRILINOS_LIBS)

ex_VbrMatrix_SOURCES = $(srcdir)/ex_VbrMatrix.cpp
ex_VbrMatrix_CXXFLAGS = $(WEBTRILINOS_INCLUDES) 
ex_VbrMatrix_LDADD = $(WEBTRILINOS_LIBS)
