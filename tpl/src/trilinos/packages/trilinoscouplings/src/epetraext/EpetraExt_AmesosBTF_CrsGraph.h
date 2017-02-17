//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
//                                                                                                    
#ifndef EpetraExt_AMESOS_BTF_CRSGRAPH_H
#define EpetraExt_AMESOS_BTF_CRSGRAPH_H

#include <EpetraExt_Transform.h>
#include <Teuchos_RCP.hpp>
#include <vector>

class Epetra_CrsGraph;
class Epetra_Map;

namespace EpetraExt {

///
/** Block Triangular Factorization (Reordering) of Epetra_CrsGraph
 *
 * Uses Tim Davis' BTF algorithm to find a block lower or upper triangular
 * ordering form a Epetra_CrsGraph.
 */

class AmesosBTF_CrsGraph : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 public:

  ///
  /** Destructor
   */
  ~AmesosBTF_CrsGraph();

  ///
  /** Default Constructor
   */
  AmesosBTF_CrsGraph( bool upperTri = false, bool verbose = false )
    : numBlocks_(0),
      upperTri_(upperTri),
      verbose_(verbose)
    {}
    
  ///
  /** Construction of BTF ordered Epetra_CrsGraph from <tt>orig</tt> object.
   *
   * Preconditions:<ul>
   * </ul>
   *
   * Invariants:<ul>
   * </ul>
   *
   * Postconditions:<ul>
   * </ul>
   *
   */
  NewTypeRef operator()( OriginalTypeRef orig );

  std::vector<int> RowPerm() { return rowPerm_; }
  std::vector<int> ColPerm() { return colPerm_; }
  std::vector<int> BlockPtr() { return blkPtr_; }
  int NumBlocks() { return numBlocks_; }

 private:

  Teuchos::RCP<Epetra_Map> NewRowMap_;
  Teuchos::RCP<Epetra_Map> NewColMap_;
  
  Teuchos::RCP<Epetra_CrsGraph> NewGraph_;

  std::vector<int> rowPerm_, colPerm_, blkPtr_;

  int numBlocks_;

  const bool upperTri_, verbose_;
};

} //namespace EpetraExt

#endif //EpetraExt_AMESOS_BTF_CRSGRAPH_H
