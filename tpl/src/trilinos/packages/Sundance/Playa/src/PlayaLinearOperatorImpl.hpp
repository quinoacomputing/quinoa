/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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

#ifndef PLAYA_LINEAROPERATORIMPL_HPP
#define PLAYA_LINEAROPERATORIMPL_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaInverseOperatorDecl.hpp"
#include "PlayaSimpleTransposedOpDecl.hpp"
#include "PlayaBlockOperatorBaseDecl.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaOut.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#endif



using namespace Playa;
using namespace Teuchos;


template <class Scalar>
class InverseOperator;


//=======================================================================
template <class Scalar>
LinearOperator<Scalar>::LinearOperator() 
  : Handle<LinearOperatorBase<Scalar> >() {;}


//=======================================================================
template <class Scalar>
LinearOperator<Scalar>::LinearOperator(const RCP<LinearOperatorBase<Scalar> >& smartPtr) 
  : Handle<LinearOperatorBase<Scalar> >(smartPtr)  {;}




//=======================================================================
template <class Scalar> inline 
void LinearOperator<Scalar>::apply(const Vector<Scalar>& in,
  Vector<Scalar>& out) const
{
  Tabs tab(0);
  PLAYA_MSG1(this->verb(), tab << "Operator=" << this->description()
    << ",  calling apply() function");
  Tabs tab1;
  
  if (this->verb() > 2)
  {
    Tabs tab2;
    Out::os() << tab2 << "input vector = " << in << std::endl;
  }
  else if (this->verb() > 1)
  {
    Tabs tab2;
    Out::os() << tab2 << "input vector = " << in.description() << std::endl;
  }

  /* the result vector might not be initialized. If it's null,
   * create a new vector in the range space */
  if (out.ptr().get()==0)
  {
    Tabs tab2;
    PLAYA_MSG3(this->verb(), tab2 << "allocating output vector");
    out = this->range().createMember();
  }
  else
  {
    Tabs tab2;
    PLAYA_MSG3(this->verb(), tab2 << "using preallocated output vector");
  }

  this->ptr()->apply(Teuchos::NO_TRANS, in, out);

  if (this->verb() > 2)
  {
    Tabs tab2;
    Out::os() << tab2 << "output vector = " << out << std::endl;
  }
  else if (this->verb() > 1)
  {
    Tabs tab2;
    Out::os() << tab2 << "output vector = " << out.description() << std::endl;
  }

  PLAYA_MSG1(this->verb(), tab << "Operator=" << this->description()
    << ",  done with apply() function");
  
}




//=======================================================================
template <class Scalar> inline 
void LinearOperator<Scalar>::applyTranspose(const Vector<Scalar>& in,
  Vector<Scalar>& out) const
{
  Tabs tab(0);
  PLAYA_MSG1(this->verb(), tab << "Operator=" << this->description()
    << ",  calling applyTranspose() function");
  Tabs tab1;
  
  if (this->verb() > 2)
  {
    Tabs tab2;
    Out::os() << tab2 << "input vector = " << in << std::endl;
  }
  else if (this->verb() > 1)
  {
    Tabs tab2;
    Out::os() << tab2 << "input vector = " << in.description() << std::endl;
  }


  /* the result vector might not be initialized. If it's null,
   * create a new vector in the range space */
  if (out.ptr().get()==0)
  {
    Tabs tab2;
    PLAYA_MSG3(this->verb(), tab2 << "allocating output vector");
    out = this->domain().createMember();
  }
  else
  {
    Tabs tab2;
    PLAYA_MSG3(this->verb(), tab2 << "using preallocated output vector");
  }

  this->ptr()->apply(Teuchos::TRANS, in, out);

  
  
  if (this->verb() > 2)
  {
    Tabs tab2;
    Out::os() << tab2 << "output vector = " << out << std::endl;
  }
  else if (this->verb() > 1)
  {
    Tabs tab2;
    Out::os() << tab2 << "output vector = " << out.description() << std::endl;
  }

  PLAYA_MSG1(this->verb(), tab << "Operator=" << this->description()
    << ",  done with applyTranpose() function");
  
}


//=======================================================================
template <class Scalar>
RCP<Time>& LinearOperator<Scalar>::opTimer()
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("Low-level vector operations");
  return rtn;
}

//=======================================================================
template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::transpose() const
{
  LinearOperator<Scalar> op = transposedOperator(*this);
  return op;
}





//=======================================================================
template <class Scalar>
RCP<LoadableMatrix<Scalar> > LinearOperator<Scalar>::matrix()
{
  RCP<LoadableMatrix<Scalar> > rtn 
    = rcp_dynamic_cast<LoadableMatrix<Scalar> >(this->ptr());
  return rtn;
}

//=======================================================================
template <class Scalar>
void LinearOperator<Scalar>::getRow(const int& row, 
  Teuchos::Array<int>& indices, 
  Teuchos::Array<Scalar>& values) const
{
  const RowAccessibleOp<Scalar>* val = 
    dynamic_cast<const RowAccessibleOp<Scalar>* >(this->ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(val == 0, std::runtime_error, 
    "Operator not row accessible; getRow() not defined.");
  val->getRow(row, indices, values);
}

//=============================================================================
template <class Scalar>
int LinearOperator<Scalar>::numBlockRows() const
{
  const BlockOperatorBase<Scalar>* b = dynamic_cast<const BlockOperatorBase<Scalar>* >(this->ptr().get());
  if (b==0) return 1;
  return b->numBlockRows(); 
}

//=============================================================================
template <class Scalar>
int LinearOperator<Scalar>::numBlockCols() const
{
  const BlockOperatorBase<Scalar>* b = dynamic_cast<const BlockOperatorBase<Scalar>* >(this->ptr().get());
  if (b==0) return 1;
  return b->numBlockCols(); 
}


//=============================================================================
template <class Scalar>
const VectorSpace<Scalar> 
LinearOperator<Scalar>::range() const
{return this->ptr()->range();}
  

//=============================================================================
template <class Scalar>
void LinearOperator<Scalar>::setBlock(int i, int j, 
  const LinearOperator<Scalar>& sub) 
{
  SetableBlockOperatorBase<Scalar>* b = 
    dynamic_cast<SetableBlockOperatorBase<Scalar>* >(this->ptr().get());
  
  TEUCHOS_TEST_FOR_EXCEPTION(b == 0, std::runtime_error, 
    "Can't call setBlock since operator not SetableBlockOperatorBase");

  b->setBlock(i, j, sub);
} 



//=============================================================================
template <class Scalar>
const  VectorSpace<Scalar> 
LinearOperator<Scalar>::domain() const 
{return this->ptr()->domain();}



//=============================================================================
template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::getBlock(const int &i, 
  const int &j) const 
{
  const BlockOperatorBase<Scalar>* b = 
    dynamic_cast<const BlockOperatorBase<Scalar>* >(this->ptr().get());
  
  if (b==0)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(i != 0 || j != 0, std::runtime_error, 
      "nonzero block index (" << i << "," << j << ") into "
      "non-block operator");
    return *this;
  }
  return b->getBlock(i, j);
}


//=============================================================================
template <class Scalar>
LinearOperator<Scalar> LinearOperator<Scalar>::getNonconstBlock(const int &i, 
  const int &j) 
{
  BlockOperatorBase<Scalar>* b = 
    dynamic_cast<BlockOperatorBase<Scalar>* >(this->ptr().get());
  
  if (b==0)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(i != 0 || j != 0, std::runtime_error, 
      "nonzero block index (" << i << "," << j << ") into "
      "non-block operator");
    return *this;
  }
  return b->getNonconstBlock(i, j);
}

 

//=============================================================================
template <class Scalar>
void LinearOperator<Scalar>::endBlockFill() 
{
  SetableBlockOperatorBase<Scalar>* b = 
    dynamic_cast<SetableBlockOperatorBase<Scalar>* >(this->ptr().get());
  
  TEUCHOS_TEST_FOR_EXCEPTION(b == 0, std::runtime_error, 
    "Can't call endBlockFill because operator is not a SetableBlockOperator");

  
  b->endBlockFill();
} 






#endif
