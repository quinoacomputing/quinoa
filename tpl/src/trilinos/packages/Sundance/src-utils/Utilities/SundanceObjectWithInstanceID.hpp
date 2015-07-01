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

#ifndef SUNDANCE_OBJECTWITHINSTANCEID_H
#define SUNDANCE_OBJECTWITHINSTANCEID_H

#include "SundanceDefs.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
    /**
     * ObjectWithInstanceID provides a common method for the
     * generation of instance-specific ID numbers. Subclasses will inherit
     * the id() method, and instances of those subclasses will be
     * given a unique ID at construction time.
     * 
     * <h4> Design note: </h4> By templating on the derived type, 
     * we can give each derived type its own sequence of ID numbers.
     */
    template <class T>
    class ObjectWithInstanceID
    {
    public:
      /** Empty ctor will assign ID at construction time */
      ObjectWithInstanceID() : id_(nextID()) {;} 

      /** Return this object's ID number */
      int id() const {return id_;}

    private:
      /** Generate the next ID in the sequence */
      static int& nextID() {static int rtn=0; rtn++; return rtn;}

      /** */
      int id_;
    };

}


#endif  /* DOXYGEN_DEVELOPER_ONLY */   
#endif
