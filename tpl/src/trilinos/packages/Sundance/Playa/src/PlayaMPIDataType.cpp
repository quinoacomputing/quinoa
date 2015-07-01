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

#include "PlayaMPIDataType.hpp"

#ifdef HAVE_MPI
// Provide an explicit template specialization for the opaque type MPI_Datatype
// so that the instantiation of Teuchos::RCP<MPI_Datatype> objects compiles correctly in debug mode
// without relying on the implementation details of the MPI library.
#include "Teuchos_TypeNameTraits.hpp"
namespace Teuchos
{
  TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(MPI_Datatype);
} // namespace Teuchos
#endif

namespace Playa
{
using Teuchos::RCP;
using Teuchos::rcp;

MPIDataType::MPIDataType(const std::string& name)
  : name_(name)
#ifdef HAVE_MPI
  , mpiType_()
#endif
{}

#ifdef HAVE_MPI
MPIDataType::MPIDataType(const std::string& name, 
  const RCP<MPI_Datatype>& mpiType)
  : name_(name), mpiType_(mpiType){}


MPI_Datatype* MPIDataType::ptr() 
{
  TEUCHOS_TEST_FOR_EXCEPT(mpiType_.get()==0);
  return mpiType_.get();
}

const MPI_Datatype& MPIDataType::handle() const 
{
  TEUCHOS_TEST_FOR_EXCEPT(mpiType_.get()==0);
  return *(mpiType_.get());
}
#endif



std::stack<MPIDataType>& MPIDataType::typeRegistry()
{
  static std::stack<MPIDataType> rtn;

  return rtn;
}

void MPIDataType::clearTypeRegistry()
{
  while(!typeRegistry().empty())
  {
#ifdef HAVE_MPI
//    MPIDataType t = typeRegistry().top();
//    std::cerr << "clearing type " << typeRegistry().top().name() << std::endl;

    int ierr = MPI_Type_free(typeRegistry().top().ptr());

    TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::runtime_error,
      "Error code=" << ierr << " detected in MPI_Type_free()");
#endif
    typeRegistry().pop();
  }
}

void MPIDataType::registerType(const MPIDataType& dataType)
{
  typeRegistry().push(dataType);
}

MPIDataType MPIDataType::intType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI int", rcp(new MPI_Datatype(MPI_INT)));
#else
  static MPIDataType rtn("MPI int");
#endif
  return rtn;
}



MPIDataType MPIDataType::floatType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI float", rcp(new MPI_Datatype(MPI_FLOAT)));
#else
  static MPIDataType rtn("MPI float");
#endif
  return rtn;
}



MPIDataType MPIDataType::doubleType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI double", rcp(new MPI_Datatype(MPI_DOUBLE)));
#else
  static MPIDataType rtn("MPI double");
#endif
  return rtn;
}



MPIDataType MPIDataType::doubleIntPairType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI double/int pair", rcp(new MPI_Datatype(MPI_DOUBLE_INT)));
#else
  static MPIDataType rtn("MPI double/int pair");
#endif
  return rtn;
}




MPIDataType MPIDataType::charType()
{
#ifdef HAVE_MPI
  static MPIDataType rtn("MPI char", rcp(new MPI_Datatype(MPI_CHAR)));
#else
  static MPIDataType rtn("MPI char");
#endif
  return rtn;
}



	
} // namespace Playa

