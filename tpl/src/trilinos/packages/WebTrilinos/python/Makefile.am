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

EXTRA_DIST = index.html help.html \
	     setpath.py \
             ex_Amesos.py  ex_Comm.py \
	     ex_Map.py     ex_ML.py \
	     ex_AztecOO.py ex_IFPACK.py \
	     ex_Matrix.py  ex_MultiVector.py  

EXAMPLE_SUITE = ex_Amesos  ex_Comm\
                ex_Map     ex_ML \
                ex_AztecOO ex_IFPACK \
                ex_Matrix  ex_MultiVector  

SETPATH = setpath.py

all: $(EXAMPLE_SUITE)
	        if [ $(srcdir) != '.' ]; then cp $(srcdir)/$(SETPATH) . ; fi

ex%: $(srcdir)/ex%.py
	        if [ $(srcdir) != '.' ]; then cp $(srcdir)/$@.py . ; fi

CLEANFILES = $(EXAMPLE_SUITE:=.py) $(SETPATH) *.pyc
