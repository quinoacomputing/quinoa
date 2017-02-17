// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef GALERI_RANDOM_H
#define GALERI_RANDOM_H

#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Util.h"
#include <limits>

namespace Galeri {
namespace Maps {

template<typename int_type>
inline
Epetra_Map* 
TRandom(const Epetra_Comm& Comm, const int_type n)
{
  if (n <= 0 || n > std::numeric_limits<int>::max())
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::Random()",
                    "n = " + toString(n)));
  
  // This is the idea: I create the map on proc 0, then I broadcast
  // it to all procs. This is not very efficient, but saves some MPI calls.
      
  std::vector<int> part(n);
      
  if (Comm.MyPID() == 0) 
  {
    Epetra_Util Util;

    for (int_type i = 0 ; i < n ; ++i) 
    {
      unsigned int r = Util.RandomInt();	
      part[i] = r%(Comm.NumProc());
    }
  }
      
  Comm.Broadcast(&part[0], n, 0);
      
  // count the elements assigned to this proc
  int NumMyElements = 0;
  for (int_type i = 0 ; i < n; ++i) 
    if (part[i] == Comm.MyPID()) NumMyElements++;

  // get the loc2global list
  std::vector<int_type> MyGlobalElements(NumMyElements);
  int count = 0;
  for (int_type i = 0 ; i < n ; ++i)
    if (part[i] == Comm.MyPID()) 
      MyGlobalElements[count++] = i;

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  return(new Epetra_Map(n, NumMyElements, (NumMyElements ? &MyGlobalElements[0] : NULL), 0, Comm));
#else
  return(new Epetra_Map);
#endif

} // TRandom()

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map* 
Random(const Epetra_Comm& Comm, const int n) {
  return TRandom<int>(Comm, n);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map* 
Random64(const Epetra_Comm& Comm, const long long n) {
  return TRandom<long long>(Comm, n);
}
#endif

} // namespace BlockMaps
} // namespace Galeri
#endif
