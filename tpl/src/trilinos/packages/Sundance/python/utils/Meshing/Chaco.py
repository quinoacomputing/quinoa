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
# Python wrapper for Chaco partitioning 
#
############################################################################

from sets import Set
from Mesh import Mesh
import posix
import string

class Chaco :

    # create a partitioner object
    def __init__(self, filename, indexOffset) :
        self.filename_ = filename
        self.indexOffset_ = indexOffset


    # driver routine for the partitioning
    def partition(self, mesh, nProc) :

        print 'creating graph...'
        mesh.writeGraph(self.filename_)
        self.runChaco(nProc)
        return self.readElemAssignments(self.filename_)

    # Set up chaco input file and run chaco

    def runChaco(self, nProcs) :
        paramsFile = file('User_Params', 'w')
        paramsFile.write('OUTPUT_ASSIGN=true\n')
        paramsFile.write('PROMPT=false\n')
        paramsFile.write('ARCHITECTURE=1\n')
        paramsFile.write('REFINE_PARTITION=4\n')
        paramsFile.write('REFINE_MAP=true\n')
        paramsFile.write('KL_BAD_MOVES=20\n')
        paramsFile.write('KL_NTRIES_BAD=10\n')
        paramsFile.write('KL_IMBALANCE=0.02\n')
        paramsFile.write('INTERNAL_VERTICES=true\n')
        paramsFile.write('MATCH_TYPE=2\n')
        paramsFile.write('HEAVY_MATCH=true\n')
        paramsFile.write('TERM_PROP=true\n')
        paramsFile.write('COARSE_NLEVEL_KL=1\n')
        paramsFile.write('COARSEN_RATIO_MIN=0.7\n')
        paramsFile.write('CUT_TO_HOP_COST=1.0\n')
        paramsFile.write('RANDOM_SEED=12345\n')
        paramsFile.flush()
        inputFile = file('chacoInput', 'w')
        inputFile.write('%s.graph\n' % self.filename_)
        inputFile.write('%s.assign\n' % self.filename_)
        inputFile.write('1\n')    # select partitioning algorithm. 1=multilevel
        inputFile.write('100\n')  # num vertices in coarsest graph
        inputFile.write('%d\n' % nProcs) # number of partitions
        inputFile.write('1\n') # apply bisection
        inputFile.write('n\n') # answer 'n' to query about running another problem
        inputFile.flush()
        inputFile.close()
        print 'calling Chaco'
        posix.system('chaco < chacoInput')

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
            partFile.write('%d %d\n' % (i+self.indexOffset_, assignments[i]+self.indexOffset_))
