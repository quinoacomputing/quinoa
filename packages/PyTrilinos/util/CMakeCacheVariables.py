#! /usr/bin/env python

# @HEADER
# ************************************************************************
# 
#              PyTrilinos: Python Interface to Trilinos
#                 Copyright (2004) Sandia Corporation
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
CMakeVariables.py - This script can be run as a shell command or imported as a
module into a python script.  Its purpose is to read a CMake cache file and
obtain a map of CMake variables and their values.

Shell usage: CMakeCacheVariables.py [options] [filename]

options:
    -h | --help       Print this message and exit
    -c | --cmake      Output in CMake format
    -p | --python     Output as a python dictionary (default)
    -v | --version    Print the version number and exit

Python usage: import CMakeCacheVariables

Available function:

    parseCMakeCacheFile(string) -> dict
        Given the name of a CMake cache file, parse the file for CMake variable
        names and values, converting BOOL variables to python bool and skipping
        type INTERNAL.
    
"""

__version__ = 1.0
__author__  = "Bill Spotz"
__date__    = "Oct 24 2010"

import optparse
import os
import re

assignRE = re.compile(r"(.+):(.+)=(.*)$")

#############################################################################

def parseCMakeCacheFile(filename):
    """
    Open filename, read in the text and return a dictionary of CMake variable
    names and values.
    """
    absFileName = os.path.abspath(filename)
    lines = open(absFileName,"r").readlines()
    lines = [s.split("#")[0] for s in lines]
    dict = { }
    for line in lines:
        match = assignRE.match(line)
        if match:
            varName  = match.group(1).strip()
            varType  = match.group(2).strip().lower()
            varValue = match.group(3).strip()
            if varType == "bool":
                dict[varName] = varValue.lower() in ("on", "true", "yes")
            elif varType == "internal":
                pass
            else:
                dict[varName] = varValue
    return dict

#############################################################################

def main():
    """
    This is the routine that gets called if MakefileVariable.py is invoked from
    the shell.  Process any command-line options that may be given, take the
    first argument as a filename, process it to obtain a dictionary of variable
    name/value pairs, and output the results.
    """
    # Initialization
    parser = optparse.OptionParser()
    parser.set_defaults(out_type="python")
    parser.add_option("-c", "--cmake", action="store_const", const="cmake",
                      dest="out_type", help="Output results in CMake format")
    parser.add_option("-p", "--python", action="store_const", const="python",
                      dest="out_type", help="Output results in Python format")
    parser.add_option("-v", "--version", action="store_true", dest="version",
                      default=False, help="Print the version number")
    options,args = parser.parse_args()

    # Version
    if options.version:
        print "CMakeCacheVariables.py version", __version__
        exit()
    
    # Process the filename
    dict = parseCMakeCacheFile(args[0])

    # Output the variable names and values
    if options.out_type == "python":
        print dict
    elif options.out_type == "cmake":
        varNames = dict.keys()
        varNames.sort()
        for varName in varNames:
            varValue = dict[varName]
            if isinstance(varValue, bool):
                varType = "BOOL"
                if varValue:
                    varValue = "TRUE"
                else:
                    varValue = "FALSE"
            elif os.path.exists(varValue):
                if os.path.isdir(varValue):
                    varType = "PATH"
                else:
                    varType = "FILEPATH"
            else:
                varType  = "STRING"
            print "%s:%s=%s" % (varName, varType, varValue)

#############################################################################
# If called from the command line, call main()
#############################################################################

if __name__ == "__main__":
    main()
