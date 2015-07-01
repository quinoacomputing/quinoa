#@HEADER@
# ************************************************************************
# 
#                              Sundance
#                 Copyright (2005) Sandia Corporation
# 
# Copyright (year first published) Sandia Corporation.  Under the terms 
# of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
# retains certain rights in this software.
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
# Questions? Contact Kevin Long (krlong@sandia.gov), 
# Sandia National Laboratories, Livermore, California, USA
# 
# ************************************************************************
#@HEADER@


############################################################################
#
# Python wrapper for Metis partitioning 
#
############################################################################

from sets import Set
from Mesh import Mesh
import posix
import string

class Metis :

    # create a partitioner object
    def __init__(self, filename) :
        self.filename_ = filename


    # driver routine for the partitioning
    def partition(self, mesh, nProc) :

        mesh.writeGraph(self.filename_)
        self.runMetis(nProc)
        return self.readElemAssignments(self.filename_)

    # Set up metis input file and run metis

    def runMetis(self, nProcs) :
        posix.system('kmetis %s.graph %d' % (self.filename_, nProcs) )
        posix.system('cp %s.graph.part.%d %s.assign' % (self.filename_, nProcs, self.filename_) )

    def readElemAssignments(self, filename) :
        assignments = []
        f = file(self.filename_ + '.assign')

        # read lines in the assignment file
        while 1:
            line = f.readline()
            if not line : break
            if line[0]=='#': continue
            assignments.append(int(line))

        return assignments


    def writePartitionFile(self, filename, assignments, nProcs) :
        partFile = file(filename + '.part', 'w')
        partFile.write('%d %d\n' % (len(assignments), nProcs))
        for i in range(len(assignments)) :
            partFile.write('%d %d\n' % (i, assignments[i]))
