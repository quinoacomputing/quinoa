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

#include "PlayaMPIOp.hpp"

#ifdef HAVE_MPI
// Provide an explicit template specialization for the opaque type MPI_Op
// so that the instantiation of Teuchos::RCP<MPI_Op> objects compiles correctly in debug mode
// without relying on the implementation details of the MPI library.
#include "Teuchos_TypeNameTraits.hpp"
namespace Teuchos
{
  TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(MPI_Op);
} // namespace Teuchos
#endif

namespace Playa
{
using Teuchos::RCP;
using Teuchos::rcp;

MPIOp::MPIOp(const std::string& name)
  : name_(name)
#ifdef HAVE_MPI
  , mpiOp_()
#endif
{}

#ifdef HAVE_MPI
MPIOp::MPIOp(const std::string& name, 
  const RCP<MPI_Op>& mpiOp)
  : name_(name), mpiOp_(mpiOp){}


MPI_Op* MPIOp::ptr() 
{
  TEUCHOS_TEST_FOR_EXCEPT(mpiOp_.get()==0);
  return mpiOp_.get();
}

const MPI_Op& MPIOp::handle() const 
{
  TEUCHOS_TEST_FOR_EXCEPT(mpiOp_.get()==0);
  return *(mpiOp_.get());
}
#endif



std::stack<MPIOp>& MPIOp::opRegistry()
{
  static std::stack<MPIOp> rtn;

  return rtn;
}

void MPIOp::clearOpRegistry()
{
  while(!opRegistry().empty())
  {
#ifdef HAVE_MPI
    MPIOp t = opRegistry().top();

    int ierr = MPI_Op_free(t.ptr());

    TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::runtime_error,
      "Error code=" << ierr << " detected in MPI_Type_free()");
#endif
    opRegistry().pop();
  }
}

void MPIOp::registerOp(const MPIOp& opType)
{
  opRegistry().push(opType);
}


MPIOp MPIOp::sumOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI sum", rcp(new MPI_Op(MPI_SUM)));
#else
  static MPIOp rtn("MPI sum");
#endif
  return rtn;
}


MPIOp MPIOp::minOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI min", rcp(new MPI_Op(MPI_MIN)));
#else
  static MPIOp rtn("MPI min");
#endif
  return rtn;
}

MPIOp MPIOp::maxOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI max", rcp(new MPI_Op(MPI_MAX)));
#else
  static MPIOp rtn("MPI max");
#endif
  return rtn;
}

MPIOp MPIOp::minlocOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI minloc", rcp(new MPI_Op(MPI_MINLOC)));
#else
  static MPIOp rtn("MPI minloc");
#endif
  return rtn;
}

MPIOp MPIOp::maxlocOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI maxloc", rcp(new MPI_Op(MPI_MAXLOC)));
#else
  static MPIOp rtn("MPI maxloc");
#endif
  return rtn;
}

MPIOp MPIOp::productOp()
{
#ifdef HAVE_MPI
  static MPIOp rtn("MPI prod", rcp(new MPI_Op(MPI_PROD)));
#else
  static MPIOp rtn("MPI prod");
#endif
  return rtn;
}


	
} // namespace Playa

