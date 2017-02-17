#! /usr/bin/env python

# @HEADER
# ************************************************************************
# 
#              PyTrilinos: Python Interface to Trilinos
#                 Copyright (2010) Sandia Corporation
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
# 
# ************************************************************************
# @HEADER

"""
copyWithCMakeSubstitutions.py - This script is intended to be run as a script.
Its purpose is to copy a file from one location to another, and substitute
instances of '${VAR_NAME}' with the value of variable VAR_NAME, which comes from
the local build tree's CMakeCache.txt file.  If VAR_NAME is not a variable name
within the cache file, then substitute the null string.
"""

__version__ = "1.0"
__author__  = "Bill Spotz"
__date__    = "Nov 4 2010"

# Import python modules
import MakefileVariables
import CMakeCacheVariables
import optparse
import os
import re

#############################################################################

def main():
    """
    Process any command line arguments.  Query the build system to find the
    top-level binary directory, to find the CMakeCache.txt file.  Query the
    CMakeCache.txt file to get a dictionary of all the cache variables.  Read in
    the input file, make cache variable substitutions, and write the output
    file.
    """

    # Parse the command line options
    parser = optparse.OptionParser()
    parser.add_option("-v", "--version", action="store_true", dest="version",
                      default=False, help="Print the version number and exit")
    options,args = parser.parse_args()

    # Version
    if options.version:
        print "copyWithCMakeSubstitutions.py version", __version__
        exit()

    # Arguments
    if len(args) != 2:
        print "usage: copyWithCMakeSubstitutions.py infile outfile"
        exit(-1)
    infilename  = args[0]
    outfilename = args[1]
    if infilename == outfilename:
        print "outfile must be distinct from infile"
        exit(-2)

    # Obtain the dictionary of CMake cache variables
    makefileVariables   = MakefileVariables.processMakefile("Makefile")
    cmakeBinaryDir      = makefileVariables["CMAKE_BINARY_DIR"]
    cmakeCacheFile      = os.path.join(cmakeBinaryDir, "CMakeCache.txt")
    cmakeCacheVariables = CMakeCacheVariables.parseCMakeCacheFile(cmakeCacheFile)

    # Open the files
    infile  = open(infilename,  "r")
    outfile = open(outfilename, "w")

    # Read in the input source file, line by line, make substitutions, and write
    # the output target file
    substitution = re.compile(r"\$\{.*\}")
    line = infile.read()
    while line:
        match = re.search(substitution, line)
        while match:
            name  = line[match.start()+2:match.end()-1]
            line  = line.replace(match.group(), str(cmakeCacheVariables.get(name,"")))
            match = re.search(substitution, line)
        outfile.write(line)
        line = infile.read()

    # Close the files and transfer the protection mode
    infile.close()
    outfile.close()
    mode = os.stat(infilename)[0]
    os.chmod(outfilename, mode)

#############################################################################
# If called from the command line, call main()
#############################################################################

if __name__ == "__main__":
    main()
