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
# Python mesh class
#
############################################################################

from sets import Set

class Mesh :
    def __init__(self) :
        self.pts_ = []
        self.nodeAttrs_ = []
        self.elemAttrs_ = []
        self.marks_ = []
        self.elemVerts_ = []
        self.vertToElemMap_ = {}
        self.dim_ = 0
        self.numElemAttr_ = 0
        self.numNodeAttr_ = 0
        self.numMark_ = 0
        self.sides_ = {}

    def setDimension(self, d) :
        self.dim_ = d

    def setNumAttributes(self, numAttr) :
        self.numAttr_ = numAttr

    def setNumElemAttributes(self, numAttr) :
        self.numElemAttr_ = numAttr
        
    def setNumMarkers(self, numMark) :
        self.numMark_ = numMark
        
    def numPts(self) :
        return len(self.pts_)

    def numElems(self) :
        return len(self.elemVerts_)

    def numNodeAttributes(self) :
        return self.numNodeAttr_

    def numElemAttributes(self) :
        return self.numElemAttr_

    def numMarkers(self) :
        return self.numMark_

    def dim(self) :
        return self.dim_;

    def addPoint(self, pt) :
        self.pts_.append(pt)

    def addMarker(self, mark) :
        self.marks_.append(mark)

    def addNodeAttrs(self, attrs) :
        self.nodeAttrs_.append(attrs)

    def addElemAttrs(self, attrs) :
        self.elemAttrs_.append(attrs)

    def addElem(self, elemVerts) :
        lid = len(self.elemVerts_)
        self.elemVerts_.append(elemVerts)

        for v in elemVerts :
            if v not in self.vertToElemMap_ :
                self.vertToElemMap_[v] = Set()
            self.vertToElemMap_[v].add(lid)

    def addSide(self, elemGID, facet, label) :
        if elemGID in self.sides_ :
            self.sides_[elemGID].append((facet, label))
        else:
            self.sides_[elemGID] = [(facet, label)]

    def getElemSides(self, elemGID) :
        return self.sides_[elemGID]

    def elemNumSides(self, elemGID) :
        if elemGID in self.sides_:
            return len(self.sides_[elemGID])
        else:
            return 0

    def hasSides(self) :
        return len(self.sides_) > 0

    def sides(self) :
        return self.sides_


    def getNodeData(self, index) :
        ptStr = self.ptToStr(self.pts_[index])
        
        if self.numNodeAttr_==0 :
            attrStr = '';
        else:
            attrStr = self.attrsToStr(self.nodeAttrs_[index])
        
        if self.numMark_==0 :
            markStr = '';
        else:
            markStr = self.attrsToStr(self.marks_[index])

        if self.numNodeAttr_==0 :
            if self.numMark_==0 :
                return '%s' % ptStr
            else:
                return '%s %s' % (ptStr, markStr)
        else:
            if self.numMark_==0 :
                return '%s %s' % (ptStr, attrStr)
            else:
                return '%s %s %s' % (ptStr, attrStr, markStr)


    def getElemData(self, index, nodeGIDToLIDMap) :
        vertStr = self.vertsToStr(self.elemVerts_[index], nodeGIDToLIDMap)
        
        if self.numElemAttr_==0 :
            return '%s' % vertStr
        else:
            attrStr = self.elemAttrs_[index]
            return '%s %s' % (vertStr, attrStr)

    def ptToStr(self, pt):
        if self.dim_==2:
            return '%g %g' % (pt[0], pt[1])
        else:
            return '%g %g %g' % (pt[0], pt[1], pt[2])

    def vertsToStr(self, verts, nodeGIDToLIDMap):
        rtn = '';
        for v in verts :
            vs = '%d ' % nodeGIDToLIDMap[v]
            rtn = rtn + vs
        return rtn

    def attrsToStr(self, attrs) : 
        rtn = '';
        for a in attrs :
            s = '%d ' % a
            rtn = rtn + s
        return rtn

    def findNeighbors(self) :
        neighbors = []
        nEdges = 0
        for i in range(self.numElems()) :
            allNeighbors = Set()
            for v in self.elemVerts_[i] :
                allNeighbors = Set.union(allNeighbors, self.vertToElemMap_[v])
            # get rid of self-references
            allNeighbors.discard(i)
            fullNeighbors = []
            for j in allNeighbors :

                numCommonNodes = Set.intersection(self.elemVerts_[i],
                                                  self.elemVerts_[j])
                if len(numCommonNodes) == self.dim_ :
                    fullNeighbors.append(j)
                         
            nEdges = nEdges + len(fullNeighbors)
            neighbors.append(fullNeighbors)

        nEdges = nEdges/2

        return (neighbors, nEdges)

    def writeGraph(self, filename) :
         print 'determining neighbors...'
         (neighbors, nEdges) = self.findNeighbors()
         graphFile = file(filename + '.graph', 'w')
         print 'writing graph file...'
         graphFile.write('%d %d\n' % (self.numElems(), nEdges))

         for i in range(self.numElems()) :
             line = ''
             for j in neighbors[i] :
                 line = line +  '%d ' % (j+1) 
             graphFile.write(line + '\n');


    def remapEntities(self, assignments, nProc) :

        procBuckets = []

        for p in range(nProc) :
            procBuckets.append([])

        k = 0
        for a in assignments :
            procBuckets[a].append(k)
            k = k + 1

        g = 0
        entityMap = range(len(assignments))
        for p in range(nProc) :
            for i in procBuckets[p] :
                entityMap[i] = g
                g = g + 1

        return entityMap

    def getNodeAssignments(self, nProc, elemAssignments) :
        nodeAssignments = []
        nodeOwners = []
        nodesPerProc = [0 for i in range(nProc)]
        for node in self.vertToElemMap_.keys() :
            owner = max(self.vertToElemMap_[node])
            nodeOwners.append(owner)
            nodeAssignments.append(elemAssignments[owner])
            nodesPerProc[elemAssignments[owner]] = nodesPerProc[elemAssignments[owner]]+1
    
        return (nodeAssignments, nodeOwners, nodesPerProc)


    def getElemsPerProc(self, nProc, elemAssignments) :
        elemsPerProc = [0 for i in range(nProc)]
        for e in elemAssignments:
            elemsPerProc[e] = elemsPerProc[e]+1
        return elemsPerProc


    def getOffProcData(self, p, elemAssignments, nodeAssignments) :

        # first pass: find off-proc nodes required by on-proc elems
        offProcNodes = Set()
        for e in range(len(self.elemVerts_)) :
            if elemAssignments[e] != p :
                continue
            for v in self.elemVerts_[e] :
                if nodeAssignments[v] == p :
                    continue
                offProcNodes.add(v)

        # second pass: find off-proc elems required by the nodes
        # obtained in the previous step
        offProcElems = Set()
        for v in offProcNodes :
            for e in self.vertToElemMap_[v] :
                if elemAssignments[e] == p :
                    continue
                offProcElems.add(e)

        # third pass: find the additional nodes required by the off-proc
        # elems found in the previous step

        for e in offProcElems :
            for v in self.elemVerts_[e] :
                if nodeAssignments[v] == p :
                    continue
                offProcNodes.add(v)

        return (offProcNodes, offProcElems)
    
            
    
    def getNodeParInfo(self, p, nProc, nodeAssignments, offProcNodes) :

        nodeGID = []
        nodeOwners = []

        for n in range(len(nodeAssignments)) :
            if nodeAssignments[n]==p :
                nodeGID.append(n)
                nodeOwners.append(p)

        for n in offProcNodes :
            nodeGID.append(n)
            nodeOwners.append(nodeAssignments[n])    

        return (nodeGID, nodeOwners)
    
    def getElemParInfo(self, p, nProc, elemAssignments, offProcElems) :

        elemGID = []
        elemOwners = []

        for n in range(len(elemAssignments)) :
            if elemAssignments[n]==p :
                elemGID.append(n)
                elemOwners.append(p)

        for n in offProcElems :
            elemGID.append(n)
            elemOwners.append(elemAssignments[n])

        return (elemGID, elemOwners)
    
                
        
                
            

        
        

    
        
        
