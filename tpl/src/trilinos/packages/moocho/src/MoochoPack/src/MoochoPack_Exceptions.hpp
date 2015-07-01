// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H
#define REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H

#include "MoochoPack_Types.hpp"
#include "ConstrainedOptPack_QPSolverStats.hpp"

namespace MoochoPack {

/** \defgroup MoochoPack_grp Standard exceptions for MoochoPack */
//@{

// Thrown if the constraints are infeasible
class InfeasibleConstraints : public std::logic_error
{public: InfeasibleConstraints(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if a line search failure occurs.
class LineSearchFailure : public std::runtime_error
{public: LineSearchFailure(const std::string& what_arg) : std::runtime_error(what_arg){}};

/// Thrown if a runtime test failed.
class TestFailed : public std::runtime_error
{public: TestFailed(const std::string& what_arg) : std::runtime_error(what_arg){}};

/// Thrown if a the QP failed and was not corrected
class QPFailure : public std::runtime_error
{
public:
  QPFailure(const std::string& what_arg
        , const ConstrainedOptPack::QPSolverStats& _qp_stats)
    : std::runtime_error(what_arg)
    , qp_stats(_qp_stats)
    {}
  ConstrainedOptPack::QPSolverStats qp_stats;
};

//@}

}	// end namespace MoochoPack 

#endif // REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H
