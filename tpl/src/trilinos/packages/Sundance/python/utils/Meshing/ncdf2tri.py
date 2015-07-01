#!/usr/bin/python

# convert an NCDF file to Triangle format
# 

import sys
import posix
from sets import Set

args = Set()

def offsetInt(str) :
    return int(str)-1

def requireKeyword(f, name) :
    while 1:
        line = f.readline()
        if line[0]=="#" :
            continue

        tokens = line.split()
        if len(tokens)==0 : continue
        if tokens[0]==name: break

        raise Exception('expected keyword %s, found line %s' % (name, line))

def waitForKeyword(f, name) :
    while 1:
        line = f.readline()
        if line[0]=="#" :
            continue


        tokens = line.split()
        if len(tokens)==0 : continue
        if tokens[0]==name: break
    

def getDataBlock(f, tokens, keyword, conversionFunc) :
    rtn = []
    doneWithData = False
    doneWithBlock = False
    for t in tokens :
        if t=='=' or t==keyword : continue
        if t==';' :
            doneWithBlock=True
            break
        rtn.append(conversionFunc(t.strip(',')))
    while not doneWithBlock :
        line = f.readline()
        if line[0]=="#" :
            continue
        if line=='' :
            doneWithData = True
            break
        tokens = line.split()
        for t in tokens :
            if t=='=' or t==keyword : continue
            if t==';' :
                doneWithBlock=True
                break
            rtn.append(conversionFunc(t.strip(',')))
    return (rtn, doneWithData)

for x in sys.argv :
    args.add(x)

if (len(args) < 2) | ('-h' in args) | ('--h' in args) | ('--help' in args) :
  print 'usage: ncdf2tri <filename>'
  print 'do not include .ncdf suffix on filename'
  print 'output will be written to <filename>.ele, <filename>.node, etc'
  sys.exit(0)


filename = sys.argv[1]


f = file(filename + '.ncdf')    

# read netcdf block
requireKeyword(f, 'netcdf')

# read dimensions block
requireKeyword(f, 'dimensions:')

dimension = 0
nNodes = 0
nElem =  0
nElemBlocks = 0
nSideSets = 0
nNodeSets = 0
blockSizes = []
nodesPerElem = []
sideSetSizes = []
nodeSetSizes = []

while 1:
    line = f.readline()
    if line[0]=="#" :
        continue
    tokens = line.split()
    
    if tokens[0]=='variables:' : break
    
    if tokens[1] != '=' :
        raise Exception('expected [=] as second token, line was %s' % line)
    keyword = tokens[0]
    
    if keyword=='num_dim' :
        dimension = int(tokens[2])
    elif keyword=='num_nodes' :
        nNodes = int(tokens[2])
    elif keyword=='num_elem':
        nElem = int(tokens[2])
    elif keyword=='num_el_blk':
        nElemBlocks = int(tokens[2])
        blockSizes = [0 for i in range(nElemBlocks)]
        nodesPerElem = [0 for i in range(nElemBlocks)]
    elif keyword=='num_side_sets' :
        nSideSets = int(tokens[2])
        sideSetSizes = [0 for i in range(nSideSets)]
    elif keyword=='num_node_sets' :
        nNodeSets = int(tokens[2])
        nodeSetSizes = [0 for i in range(nNodeSets)]
    else :
        for i in range(nElemBlocks) :
            s = 'num_el_in_blk%d' % (i+1)
            if keyword==s :
                blockSizes[i] = int(tokens[2])
                s = 'num_nod_per_el%d' % (i+1)
                if keyword==s :
                    nodesPerElem[i] = int(tokens[2])
                    
        for i in range(nSideSets) :
            s = 'num_side_ss%d' % (i+1)
            if keyword==s :
                sideSetSizes[i] = int(tokens[2])
                    
        for i in range(nNodeSets) :
            s = 'num_node_ns%d' % (i+1)
            if keyword==s :
                nodeSetSizes[i] = int(tokens[2])
            
    # go forward until we find the data keyword
waitForKeyword(f, 'data:')

connect = [ [] for b in range(nElemBlocks) ]
sideSetFacets = [ [] for s in range(nSideSets) ]
sideSetElems = [ [] for s in range(nSideSets) ]
nodeSetNodes = [ [] for s in range(nNodeSets) ]
doneWithData = False

while not doneWithData :
    line = f.readline()
    if line=='' :
        doneWithData = True
        break
    if line[0]=="#" :
        continue
        
    tokens = line.split()
    if len(tokens)==0 : continue

    # read coordinates
    if tokens[0]=='coord' :
        (coords, doneWithData) = getDataBlock(f, tokens, 'coord', float)
        continue

    # read element connectivity
    for b in range(nElemBlocks) :
        keyword = 'connect%d' % (b+1)
        if tokens[0]==keyword :
            (connect[b], doneWithData) = getDataBlock(f, tokens, keyword,
                                                      offsetInt)

    # read side set element numbers
    for s in range(nSideSets) :
        keyword = 'elem_ss%d' % (s+1)
        if tokens[0]==keyword :
            (sideSetElems[s], doneWithData) = getDataBlock(f, tokens,
                                                           keyword,
                                                           offsetInt)

    # read side set facet numbers
    for s in range(nSideSets) :
        keyword = 'side_ss%d' % (s+1)
        if tokens[0]==keyword :
            (sideSetFacets[s], doneWithData) = getDataBlock(f, tokens,
                                                            keyword,
                                                            offsetInt)

    # read node set node numbers
    for s in range(nNodeSets) :
        keyword = 'node_ns%d' % (s+1)
        if tokens[0]==keyword :
            (nodeSetNodes[s], doneWithData) = getDataBlock(f, tokens,
                                                           keyword,
                                                           offsetInt)





# write the node file
        
nodefile = file('%s.node' % filename, 'w')

nodefile.write('# created from %s.ncdf\n' % filename)
nodefile.write('%d %d 0 0\n' % (nNodes, dimension))

if dimension==2 :
    for n in range(nNodes):
        nodefile.write('%d %22.16g %22.16g\n' % (n, coords[n], coords[nNodes+n]))
else :
    for n in range(nNodes):
        nodefile.write('%d %22.16g %22.16g %22.16g\n' % (n, coords[n],
                                         coords[nNodes+n], coords[2*nNodes+n]))



# write the element file

    
elemfile = file('%s.ele' % filename, 'w')

elemfile.write('# created from %s.ncdf\n' % filename)
elemfile.write('%d %d 1\n' % (nElem, dimension+1))

count = 0
for b in range(nElemBlocks):
    if dimension==2 :
        for n in range(blockSizes[b]):
            elemfile.write('%d %d %d %d %d\n' % (count, connect[b][3*n],
                                                 connect[b][3*n+1],
                                                 connect[b][3*n+2], b+1))
            count = count+1
    else :
        for n in range(blockSizes[b]):
            elemfile.write('%d %d %d %d %d %d\n' % (count, connect[b][4*n],
                                                    connect[b][4*n+1],
                                                    connect[b][4*n+2],
                                                    connect[b][4*n+3], b+1))
            count = count+1



# write the sideset file

numSides = 0
for i in range(nSideSets) :
    numSides = numSides + sideSetSizes[i]


sidefile = file('%s.side' % filename, 'w')

sidefile.write('# created from %s.ncdf\n' % filename)
sidefile.write('%d\n' % numSides)

count = 0
for ss in range(nSideSets) :
    for s in range(sideSetSizes[ss]) :
        sidefile.write('%d %d %d %d\n' % (count, sideSetElems[ss][s],
                                          sideSetFacets[ss][s], ss+1))
        count = count+1



