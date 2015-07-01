#!/usr/bin/env python

#######################################################################################
#
# Partitions a Triangle file
#
#######################################################################################

import Chaco
from Chaco import *

import TriangleReader
from TriangleReader import *

import meshUtils
from meshUtils import *

import sys
import posix
from sets import Set

args = Set()
argmap = {}

for x in sys.argv :
    args.add(x)
    z = x.split('=')
    if len(z)==2 : argmap[z[0]]=z[1]


if (len(args) < 2) | ('-h' in args) | ('--h' in args) | ('--help' in args) \
       | ('--i' not in argmap) | ('--np' not in argmap) :
    print 'usage: Partition.py --np=<num processors> --i=<filename>'
    print 'do not include .ele or .node suffices on filename'
    sys.exit(0)

  


filename = argmap['--i']
np = int(argmap['--np'])

offset = 0

r = TriangleReader(filename, offset)

mesh = r.getMesh()



c = Chaco(filename, offset)

elemAssignments = c.partition(mesh, np)


(nodeAssignments, nodeOwners, nodesPerProc) = mesh.getNodeAssignments(np, elemAssignments)

elemsPerProc = mesh.getElemsPerProc(np, elemAssignments)



for procID in range(np) :
    (offProcNodes, offProcElems) = mesh.getOffProcData(procID, elemAssignments,
                                                       nodeAssignments)
    nodeGIDToLIDMap = writeNodes(filename, mesh, procID, np, nodeAssignments,
                                 nodesPerProc[procID], offProcNodes)
    writeElems(filename, mesh, procID, np, elemAssignments,
               elemsPerProc[procID], offProcElems, nodeGIDToLIDMap)

    (nodeGID, nodeOwners) = mesh.getNodeParInfo(procID, np, nodeAssignments, offProcNodes)
    (elemGID, elemOwners) = mesh.getElemParInfo(procID, np, elemAssignments, offProcElems)
    writeParFile(filename, procID, np, nodeGID, nodeOwners, elemGID, elemOwners)
    if mesh.hasSides() :
        writeSides(filename, mesh, procID, np, elemAssignments,
                   elemsPerProc[procID], offProcElems, nodeGIDToLIDMap)
    
