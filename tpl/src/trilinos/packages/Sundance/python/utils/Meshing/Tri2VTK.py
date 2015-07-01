#! /usr/bin/env python


import PySundance
import math
from PySundance import *


#############################################################################
#
# Read a Triangle file, then write it to VTK
#
#############################################################################

def main():

    # Read the mesh from a Triangle file
    filename = 'hole-0.05'
    mesher  = TriangleMeshReader(filename)
    mesh = mesher.getMesh()

    if 0:
        print 'Checking mesh consistency'
        check = mesh.checkConsistency('meshCheck')
        if check==0 :
            print 'INCONSISTENT MESH'
    
    vecType = EpetraVectorType()
    basis = Lagrange(0)
    discSpace = DiscreteSpace(mesh, basis, vecType)

    p = float(getRank()) / float(getNProc())
    f = DiscreteFunction(discSpace, p)

    writer = VTKWriter(filename)
    writer.addMesh(mesh)
    writer.addField('procID', f)
    writer.write()

    # all done!

    

if __name__ == "__main__":
    main()

    







    
    
