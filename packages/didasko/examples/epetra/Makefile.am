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


EXAMPLE = ex1  ex2  ex3  ex4  ex5  ex6  ex7  ex8  ex9  ex10 \
	  ex11 ex12 ex13 ex14 ex15 ex16 ex17 ex18 ex19 ex20 \
	  ex21 ex22 ex23 ex24 ex25
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

ex4_SOURCES = $(srcdir)/ex4.cpp
ex4_CXXFLAGS = $(EXPORT_INCLUDES) 
ex4_LDADD = $(EXPORT_LIBS)

ex5_SOURCES = $(srcdir)/ex5.cpp
ex5_CXXFLAGS = $(EXPORT_INCLUDES) 
ex5_LDADD = $(EXPORT_LIBS)

ex6_SOURCES = $(srcdir)/ex6.cpp
ex6_CXXFLAGS = $(EXPORT_INCLUDES) 
ex6_LDADD = $(EXPORT_LIBS)

ex7_SOURCES = $(srcdir)/ex7.cpp
ex7_CXXFLAGS = $(EXPORT_INCLUDES) 
ex7_LDADD = $(EXPORT_LIBS)

ex8_SOURCES = $(srcdir)/ex8.cpp
ex8_CXXFLAGS = $(EXPORT_INCLUDES) 
ex8_LDADD = $(EXPORT_LIBS)

ex9_SOURCES = $(srcdir)/ex9.cpp
ex9_CXXFLAGS = $(EXPORT_INCLUDES) 
ex9_LDADD = $(EXPORT_LIBS)

ex10_SOURCES = $(srcdir)/ex10.cpp
ex10_CXXFLAGS = $(EXPORT_INCLUDES) 
ex10_LDADD = $(EXPORT_LIBS)

ex11_SOURCES = $(srcdir)/ex11.cpp
ex11_CXXFLAGS = $(EXPORT_INCLUDES) 
ex11_LDADD = $(EXPORT_LIBS)

ex12_SOURCES = $(srcdir)/ex12.cpp
ex12_CXXFLAGS = $(EXPORT_INCLUDES) 
ex12_LDADD = $(EXPORT_LIBS)

ex13_SOURCES = $(srcdir)/ex13.cpp
ex13_CXXFLAGS = $(EXPORT_INCLUDES) 
ex13_LDADD = $(EXPORT_LIBS)

ex14_SOURCES = $(srcdir)/ex14.cpp
ex14_CXXFLAGS = $(EXPORT_INCLUDES) 
ex14_LDADD = $(EXPORT_LIBS)

ex15_SOURCES = $(srcdir)/ex15.cpp
ex15_CXXFLAGS = $(EXPORT_INCLUDES) 
ex15_LDADD = $(EXPORT_LIBS)

ex16_SOURCES = $(srcdir)/ex16.cpp
ex16_CXXFLAGS = $(EXPORT_INCLUDES) 
ex16_LDADD = $(EXPORT_LIBS)

ex17_SOURCES = $(srcdir)/ex17.cpp
ex17_CXXFLAGS = $(EXPORT_INCLUDES) 
ex17_LDADD = $(EXPORT_LIBS)

ex18_SOURCES = $(srcdir)/ex18.cpp
ex18_CXXFLAGS = $(EXPORT_INCLUDES) 
ex18_LDADD = $(EXPORT_LIBS)

ex19_SOURCES = $(srcdir)/ex19.cpp
ex19_CXXFLAGS = $(EXPORT_INCLUDES) 
ex19_LDADD = $(EXPORT_LIBS)

ex20_SOURCES = $(srcdir)/ex20.cpp
ex20_CXXFLAGS = $(EXPORT_INCLUDES) 
ex20_LDADD = $(EXPORT_LIBS)

ex21_SOURCES = $(srcdir)/ex21.cpp
ex21_CXXFLAGS = $(EXPORT_INCLUDES) 
ex21_LDADD = $(EXPORT_LIBS)

ex22_SOURCES = $(srcdir)/ex22.cpp
ex22_CXXFLAGS = $(EXPORT_INCLUDES) 
ex22_LDADD = $(EXPORT_LIBS)

ex23_SOURCES = $(srcdir)/ex23.cpp
ex23_CXXFLAGS = $(EXPORT_INCLUDES) 
ex23_LDADD = $(EXPORT_LIBS)

ex24_SOURCES = $(srcdir)/ex24.cpp
ex24_CXXFLAGS = $(EXPORT_INCLUDES) 
ex24_LDADD = $(EXPORT_LIBS)

ex25_SOURCES = $(srcdir)/ex25.cpp
ex25_CXXFLAGS = $(EXPORT_INCLUDES) 
ex25_LDADD = $(EXPORT_LIBS)

include $(top_builddir)/Makefile.export.didasko

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(DIDASKO_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(DIDASKO_LIBS))
else
EXPORT_INCLUDES = $(DIDASKO_INCLUDES)
EXPORT_LIBS = $(DIDASKO_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

