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

"""Set the python search path to include the library build directory
created by the python distutils module."""

# System includes
import os
import sys
from   distutils.util import get_platform

# Obtain the current directory name
myDir,myName = os.path.split(__file__)

# Construct the the build library directory name
libDir = "lib.%s-%s" % (get_platform(), sys.version[0:3])

# Get the full path to the build directories
fullPath   = os.path.normpath(os.path.join(myDir, "..", "src", "build", libDir,
                                           "PyTrilinos"))
epetraPath = os.path.normpath(os.path.join(myDir, "..", "..", "epetra",
                                           "python", "src", "build", libDir,
                                           "PyTrilinos"))
epetraextPath = os.path.normpath(os.path.join(myDir, "..", "..", "epetraext",
                                           "python", "src", "build", libDir,
                                           "PyTrilinos"))
galeriPath = os.path.normpath(os.path.join(myDir, "..", "..", "galeri",
                                             "python", "src", "build", libDir,
                                             "PyTrilinos"))
aztecooPath = os.path.normpath(os.path.join(myDir, "..", "..", "aztecoo",
                                            "python", "src", "build", libDir,
                                            "PyTrilinos"))
mlPath = os.path.normpath(os.path.join(myDir, "..", "..", "ml",
                                            "python", "src", "build", libDir,
                                            "PyTrilinos"))
ifpackPath = os.path.normpath(os.path.join(myDir, "..", "..", "ifpack",
                                            "python", "src", "build", libDir,
                                            "PyTrilinos"))
amesosPath = os.path.normpath(os.path.join(myDir, "..", "..", "amesos",
                                            "python", "src", "build", libDir,
                                            "PyTrilinos"))
teuchosPath = os.path.normpath(os.path.join(myDir, "..", "..", "teuchos",
                                            "python", "src", "build", libDir,
                                            "PyTrilinos"))

# Insert the full path to the build library directory
# at the beginning of the python search path
if fullPath:
    sys.path.insert(0,fullPath)
if epetraPath:
    sys.path.insert(1,epetraPath)
if epetraextPath:
    sys.path.insert(1,epetraextPath)
if galeriPath:
    sys.path.insert(2,galeriPath)
if aztecooPath:
    sys.path.insert(3,aztecooPath)
if mlPath:
    sys.path.insert(4,mlPath)
if ifpackPath:
    sys.path.insert(5,ifpackPath)
if amesosPath:
    sys.path.insert(6,amesosPath)
if teuchosPath:
    sys.path.insert(7,teuchosPath)
