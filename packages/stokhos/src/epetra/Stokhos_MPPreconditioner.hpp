// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_MP_PRECONDITIONER_HPP
#define STOKHOS_MP_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"
#include "Stokhos_BlockDiagonalOperator.hpp"
#include "Epetra_Vector.h"

namespace Stokhos {

  /*! 
   * \brief An abstract class to represent a generic stochastic Galerkin 
   * preconditioner as an Epetra_Operator.
   */
  class MPPreconditioner : public virtual Epetra_Operator {
  public:

    //! Constructor
    MPPreconditioner() {}

    //! Destructor
    virtual ~MPPreconditioner() {}

    //! Setup preconditioner
    virtual void 
    setupPreconditioner(
      const Teuchos::RCP<Stokhos::BlockDiagonalOperator>& mp_op, 
      const Epetra_Vector& x) = 0;

  }; // class MPPreconditioner

} // namespace Stokhos

#endif // STOKHOS_MP_PRECONDITIONER_HPP
