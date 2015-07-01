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

ABSTRACT_LIN_ALG_PACK_THYRA_HPP = \
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_InnerProductThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_MatrixOpNonsingThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_MatrixOpThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_MultiVectorMutableThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_ThyraAccessors.hpp\
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_VectorMutableThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_VectorSpaceFactoryThyra.hpp\
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_VectorSpaceThyra.hpp

ABSTRACT_LIN_ALG_PACK_THYRA_CPP = \
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_InnerProductThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_VectorSpaceFactoryThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_VectorSpaceThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_VectorMutableThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_MatrixOpThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_MatrixOpNonsingThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_MultiVectorMutableThyra.cpp \
  $(srcdir)/AbstractLinAlgPack/AbstractLinAlgPack_ThyraAccessors.cpp

NLP_INTERFACE_PACK_THYRA_HPP = \
  $(srcdir)/NLPInterfacePack/NLPInterfacePack_NLPThyraModelEvaluatorBase.hpp\
  $(srcdir)/NLPInterfacePack/NLPInterfacePack_NLPDirectThyraModelEvaluator.hpp\
  $(srcdir)/NLPInterfacePack/NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp

NLP_INTERFACE_PACK_THYRA_CPP = \
  $(srcdir)/NLPInterfacePack/NLPInterfacePack_NLPThyraModelEvaluatorBase.cpp\
  $(srcdir)/NLPInterfacePack/NLPInterfacePack_NLPDirectThyraModelEvaluator.cpp\
  $(srcdir)/NLPInterfacePack/NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.cpp

MOOCHO_PACK_THYRA_HPP = \
  $(srcdir)/MoochoPack/MoochoPack_MoochoThyraSolver.hpp

MOOCHO_PACK_THYRA_CPP = \
  $(srcdir)/MoochoPack/MoochoPack_MoochoThyraSolver.cpp

#
# Specify locations of include files within Trilinos
#

include $(top_builddir)/Makefile.export.moocho

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(MOOCHO_INCLUDES))
else
EXPORT_INCLUDES = $(MOOCHO_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES) -I$(top_srcdir)/../thyra/example/operator_vector

# --------------------------------------------------------------------
# Moocho library specifications
# --------------------------------------------------------------------

MOOCHO_THYRA_LIB = libmoochothyra.a

MOOCHO_THYRA_H = \
	$(ABSTRACT_LIN_ALG_PACK_THYRA_HPP) \
	$(NLP_INTERFACE_PACK_THYRA_HPP) \
  $(MOOCHO_PACK_THYRA_HPP)

libmoochothyra_a_SOURCES = \
	$(ABSTRACT_LIN_ALG_PACK_THYRA_CPP) \
	$(NLP_INTERFACE_PACK_THYRA_CPP) \
  $(MOOCHO_PACK_THYRA_CPP)

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libmoochothyra_a_AR = $(ALTERNATE_AR)

else

libmoochothyra_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

libmoochothyra_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(MOOCHO_THYRA_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(MOOCHO_THYRA_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: Moocho_Config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS =

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
