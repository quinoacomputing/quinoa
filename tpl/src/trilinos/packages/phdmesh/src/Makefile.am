# @HEADER
# ************************************************************************
#       phdMesh : Parallel Heterogneous Dynamic unstructured Mesh
#                 Copyright (2007) Sandia Corporation
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
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package phdmesh subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

ELEMENT_H = \
	$(srcdir)/../include/element/Basic_Topologies.hpp \
	$(srcdir)/../include/element/CellTopology.hpp \
	$(srcdir)/../include/element/Declarations.hpp \
	$(srcdir)/../include/element/Dimensions.hpp \
	$(srcdir)/../include/element/Hexahedron_Topologies.hpp \
	$(srcdir)/../include/element/Pyramid_Topologies.hpp \
	$(srcdir)/../include/element/Quadrilateral_Topologies.hpp \
	$(srcdir)/../include/element/Stencils.hpp \
	$(srcdir)/../include/element/Tetrahedron_Topologies.hpp \
	$(srcdir)/../include/element/Triangle_Topologies.hpp \
	$(srcdir)/../include/element/Wedge_Topologies.hpp
 
ELEMENT = \
	$(srcdir)/element/ElementDeclarations.cpp \
	$(srcdir)/element/ElementDimensions.cpp \
	$(srcdir)/element/ElementTopologies.cpp

MESH_H = \
	$(srcdir)/../include/mesh/Assemble.hpp \
	$(srcdir)/../include/mesh/BulkData.hpp \
	$(srcdir)/../include/mesh/Comm.hpp \
	$(srcdir)/../include/mesh/Entity.hpp \
	$(srcdir)/../include/mesh/EntityComm.hpp \
	$(srcdir)/../include/mesh/Field.hpp \
	$(srcdir)/../include/mesh/FieldDeclarations.hpp \
	$(srcdir)/../include/mesh/FieldData.hpp \
	$(srcdir)/../include/mesh/FieldParallel.hpp \
	$(srcdir)/../include/mesh/FieldTraits.hpp \
	$(srcdir)/../include/mesh/Gather.hpp \
	$(srcdir)/../include/mesh/Kernel.hpp \
	$(srcdir)/../include/mesh/MetaData.hpp \
	$(srcdir)/../include/mesh/Part.hpp \
	$(srcdir)/../include/mesh/Proximity.hpp \
	$(srcdir)/../include/mesh/Relation.hpp \
	$(srcdir)/../include/mesh/Types.hpp

MESH = \
	$(srcdir)/mesh/MeshAssemble.cpp \
	$(srcdir)/mesh/MeshAura.cpp \
	$(srcdir)/mesh/MeshBulkData.cpp \
	$(srcdir)/mesh/MeshComm.cpp \
	$(srcdir)/mesh/MeshEntity.cpp \
	$(srcdir)/mesh/MeshEntityComm.cpp \
	$(srcdir)/mesh/MeshField.cpp \
	$(srcdir)/mesh/MeshFieldParallel.cpp \
	$(srcdir)/mesh/MeshFieldTraits.cpp \
	$(srcdir)/mesh/MeshKernel.cpp \
	$(srcdir)/mesh/MeshMetaData.cpp \
	$(srcdir)/mesh/MeshPart.cpp \
	$(srcdir)/mesh/MeshProximity.cpp \
	$(srcdir)/mesh/MeshRebalance.cpp \
	$(srcdir)/mesh/MeshRelation.cpp \
	$(srcdir)/mesh/MeshSharing.cpp \
	$(srcdir)/mesh/MeshTypes.cpp

UTIL_H = \
	$(srcdir)/../include/util/Array.hpp \
	$(srcdir)/../include/util/ArrayPrivate.hpp \
	$(srcdir)/../include/util/Basics.hpp \
	$(srcdir)/../include/util/CSet.hpp \
	$(srcdir)/../include/util/FixedPoolAlloc.hpp \
	$(srcdir)/../include/util/IndexList.hpp \
	$(srcdir)/../include/util/NamedValue.hpp \
	$(srcdir)/../include/util/NumericEnum.hpp \
	$(srcdir)/../include/util/OctTree.hpp \
	$(srcdir)/../include/util/OctTreeOps.hpp \
	$(srcdir)/../include/util/PairIter.hpp \
	$(srcdir)/../include/util/Parallel.hpp \
	$(srcdir)/../include/util/ParallelComm.hpp \
	$(srcdir)/../include/util/ParallelIndex.hpp \
	$(srcdir)/../include/util/ParallelInputStream.hpp \
	$(srcdir)/../include/util/ParallelReduce.hpp \
	$(srcdir)/../include/util/Setv.hpp \
	$(srcdir)/../include/util/TestDriver.hpp \
	$(srcdir)/../include/util/TPI.h \
	$(srcdir)/../include/util/TPI.hpp \
	$(srcdir)/../include/util/TypeList.hpp \
	$(srcdir)/../include/util/TypeName.hpp

UTIL = \
	$(srcdir)/util/ArrayPrivate.cpp \
	$(srcdir)/util/CSet.cpp \
	$(srcdir)/util/FixedPoolAlloc.cpp \
	$(srcdir)/util/NamedValue.cpp \
	$(srcdir)/util/NumericEnum.cpp \
	$(srcdir)/util/OctTree.cpp \
	$(srcdir)/util/OctTreeOps.cpp \
	$(srcdir)/util/ParallelComm.cpp \
	$(srcdir)/util/Parallel.cpp \
	$(srcdir)/util/ParallelIndex.cpp \
	$(srcdir)/util/ParallelInputStream.cpp \
	$(srcdir)/util/ParallelReduce.cpp \
	$(srcdir)/util/Setv.cpp \
	$(srcdir)/util/TestDriver.cpp \
	$(srcdir)/util/TPI_Concurrency.c \
	$(srcdir)/util/TPI_pthreads.c \
	$(srcdir)/util/TPI_Walltime.c \
	$(srcdir)/util/TypeName.cpp

MESHIO_H = \
	$(srcdir)/../include/mesh_io/ExoII.hpp

MESHIO = \
	$(srcdir)/mesh_io/MeshIO_ExoII.cpp

TXBLAS_H = \
	$(srcdir)/../include/txblas/CR4Matrix.hpp \
	$(srcdir)/../include/txblas/cr4_mxv.h \
	$(srcdir)/../include/txblas/reduction.h \
	$(srcdir)/../include/txblas/Reduction.hpp

TXBLAS = \
	$(srcdir)/txblas/txblas_axpby.c \
	$(srcdir)/txblas/txblas_CR4Matrix.cpp \
	$(srcdir)/txblas/txblas_cr4_mxv.c \
	$(srcdir)/txblas/txblas_reduction.c

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.phdmesh

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(PHDMESH_INCLUDES))
else
EXPORT_INCLUDES = $(PHDMESH_INCLUDES)
endif


AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# New_Package library specifications
# ------------------------------------------------------------------------
PHDMESH_LIB = libphdelement_mesh.a libphdmesh_util.a libphdmesh_mesh.a \
	libphdmesh_meshio.a libphdmesh_txblas.a


#np# replace new_package with the name of the package being autotool'ed here
PHDMESH_H = \
	$(ELEMENT_H) \
	$(UTIL_H) \
	$(MESH_H) \
	$(MESHIO_H) \
	$(TXBLAS_H)

libphdelement_mesh_a_SOURCES = $(ELEMENT)

libphdmesh_util_a_SOURCES = $(UTIL)

libphdmesh_mesh_a_SOURCES = $(MESH)

libphdmesh_meshio_a_SOURCES = $(MESHIO)

libphdmesh_txblas_a_SOURCES = $(TXBLAS)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libphdelement_mesh_a_AR = $(ALTERNATE_AR)
libphdmesh_util_a_AR = $(ALTERNATE_AR)
libphdmesh_mesh_a_AR = $(ALTERNATE_AR)
libphdmesh_meshio_a_AR = $(ALTERNATE_AR)
libphdmesh_txblas_a_AR = $(ALTERNATE_AR)

else

libphdelement_mesh_a_AR = $(AR) cru
libphdmesh_util_a_AR = $(AR) cru
libphdmesh_mesh_a_AR = $(AR) cru
libphdmesh_meshio_a_AR = $(AR) cru
libphdmesh_txblas_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libphdelement_mesh_a_LIBADD = $(XTRALDADD)
libphdmesh_util_a_LIBADD = $(XTRALDADD)
libphdmesh_mesh_a_LIBADD = $(XTRALDADD)
libphdmesh_meshio_a_LIBADD = $(XTRALDADD)
libphdmesh_txblas_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = $(PHDMESH_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(PHDMESH_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: phdmesh_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(includedir)
nodist_execinclude_HEADERS = phdmesh_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

