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

#ifndef SUNDANCE_CELLREORDERERIMPLEMBASE_H
#define SUNDANCE_CELLREORDERERIMPLEMBASE_H


#include "SundanceDefs.hpp"
#include "SundanceNoncopyable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include <typeinfo>

namespace Sundance
{
class MeshBase;

/**
 * Abstract interface for the low-level objects that 
 * implement cell reordering. 
 * 
 * <h4> Adding a new reordering algorithm </h4>
 *
 * To add a new reordering algorithm, you should create a new
 * subclass of CellReordererImplemBase. The only method you will
 * need to implement is
 * \code
 * virtual int advance(int currentLID) const 
 * \endcode
 * which should provide the maximal cell LID found after
 * the <tt>currentLID.</tt>
 * Depending on the algorithm , you may also want to override
 * the methods
 * \code
 * virtual int begin() const 
 * virtual int end() const 
 * \endcode
 * which return the index of the first cell to be processed,
 * and a past-the-end index. 
 */
class CellReordererImplemBase 
  : public ObjectWithClassVerbosity<CellReordererImplemBase>
{
public:
  /** Construct with a pointer to a mesh */
  CellReordererImplemBase(const MeshBase* mesh);
      
  /** virtual dtor */
  virtual ~CellReordererImplemBase(){;}

  /** return a descriptive std::string */
  virtual std::string typeName() const {return typeid(*this).name();}
    
  /** */
  virtual int advance(int currentLID) const = 0 ;
      
  /** */
  virtual int begin() const {return 0;}
      
  /** */
  virtual int end() const ;
protected:
  /** */
  const MeshBase* mesh() const {return mesh_;}

private:
  /** Unmanaged pointer to a mesh. The mesh will contain a smart
   * pointer to this reorderer, so to avoid closed reference
   * graphs we store a raw pointer here.*/
  const MeshBase* mesh_;
      
};

}



#endif
