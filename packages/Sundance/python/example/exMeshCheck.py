#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance


import math
from PySundance import *

def main():
  vecType = EpetraVectorType()
  npx = 1
  npy = getNProc()
  ny = 16
  nx = 16
  mesher  = PartitionedRectangleMesher(0.0, 1.0, nx, npx,
                                       0.0, 2.0, ny/npy, npy);
  mesh = mesher.getMesh();

  check = mesh.checkConsistency('meshCheck')
  if check==0 :
      print 'INCONSISTENT MESH'
      print 'mesh check FAILED'
  else :
      print 'mesh test PASSED'
