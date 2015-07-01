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

#ifndef PLAYA_MPITRAITS_H
#define PLAYA_MPITRAITS_H

/*! \file PlayaMPITraits.hpp
 * \brief Declaration of a templated traits class for binding MPI types to
 * C++ types. This is for use with the MPIComm class and is supposed to compile
 * rgeardless of whether MPI has been enabled. If you need to convert directly to 
 * MPI types (e.g., MPI_INT), please refer to Teuchos_MPIRawTraits.hpp.
*/

#include "PlayaMPIComm.hpp"

namespace Playa
{
	using std::string;

	/** \ingroup MPI 
	 * \brief Templated traits class that binds MPI types to C++ types
	 * \note Template specializations exist for datatypes: <tt>char</tt>,
		<tt>int</tt>, <tt>float</tt>, and <tt>double</tt>.
	 */
	template <class T> class MPITraits
		{
		public:
			/** \brief Return the MPI data type of the template argument */
			static MPIDataType type();
		};

	/** \ingroup MPI 
	 */
	template <> class MPITraits<int>
		{
		public:
			/** return the MPI data type of the template argument */
			static MPIDataType type() {return MPIDataType::intType();}
		};
	
	/** \ingroup MPI 
	 */
	template <> class MPITraits<float>
		{
		public:
			/** return the MPI data type of the template argument */
			static MPIDataType type() {return MPIDataType::floatType();}
		};
	
	/** \ingroup MPI 
	 */
	template <> class MPITraits<double>
		{
		public:
			/** return the MPI data type of the template argument */
			static MPIDataType type() {return MPIDataType::doubleType();}
		};
	
	/** \ingroup MPI 
	 */
	template <> class MPITraits<char>
		{
		public:
			/** return the MPI data type of the template argument */
			static MPIDataType type() {return MPIDataType::charType();}
		};

	
} // namespace Playa

#endif
