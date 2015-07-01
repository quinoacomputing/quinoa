/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */


#ifndef SUNDANCE_MESHIO_UTILS_HPP
#define SUNDANCE_MESHIO_UTILS_HPP

#include "Sundance.hpp"
#include "SundanceFileIOChacoPartitioner.hpp"
#include "SundanceCellLIDMappedFieldWrapper.hpp"

namespace Sundance
{

/** Read the node-based fields from a mesher */
Expr readNodalFields(const MeshSource& mesher, const Mesh& mesh,
  const VectorType<double>& vecType, int verb=0);


/** Read a 2D field stored in simple ASCII format
 * 
 * nx ny                 # num x pts, num y pts
 * u1_1 u2_1 ... uN_1    # data at node 1
 * u1_2 u2_2 ... uN_2    # data at node 2
 * 
*/
Expr readSerialGridField(const std::string& gridFile, 
  double ax, double bx,
  double ay, double by,
  const VectorType<double>& vecType,
  const MeshType& meshType,
  Mesh& mesh);


/** This function reads in an exodus file, discretizes a function on 
 * it, and computes some functional on it. It then writes the field, 
 * reads it back, and computes the same functional on the new copy. 
 * The return value is the difference between the two functionals
 * which should be zero if all is running correctly. */
double readbackTester(const std::string& infile, const MPIComm& comm, int verb=0) ;


/** Partition a mesh and write the parts to exodus files */
void serialPartition(
  const RCP<SerialPartitionerBase>& part,
  int numProc,
  const MeshSource& mesher, 
  const std::string& outfile);
}

#endif
