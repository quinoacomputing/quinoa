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

#include <string>
#include <sstream>
 
#include "IterationPack_CastIQMember.hpp"
#include "Teuchos_Assert.hpp"

namespace IterationPack {

CastIQMemberBase::CastIQMemberBase( const std::string iq_name )
  :  iq_name_(iq_name), iq_id_(NOT_SET_YET)
{}

const std::string&
CastIQMemberBase::iq_name() const
{
  return iq_name_;
}

bool
CastIQMemberBase::exists_in( const AlgorithmState& s ) const
{
  cache_iq_id(s);
  return iq_id_ != AlgorithmState::DOES_NOT_EXIST;
}

void CastIQMemberBase::cache_iq_id( const AlgorithmState& s ) const
{
  if( iq_id_ == NOT_SET_YET ) {
    iq_id_ = s.get_iter_quant_id( iq_name_ );
  }
}

void CastIQMemberBase::throw_cast_error( const AlgorithmState::iq_id_type iq_id, const std::string& iqa_name ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    iq_id == AlgorithmState::DOES_NOT_EXIST, AlgorithmState::DoesNotExist
    ,"CastIQMember<T>::operator()(...) : Error, the iteration quantity \""
    << iq_name_ << "\" does not exist in this state object." );
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, IterationPack::InvalidTypeCastException
    ,"CastIQMember<T>::operator()(state) : Error, the iteration quantity \""
    << iq_name_ << "\" exists but it is not of the type IterQuantityAccess<"
    << iqa_name << ">" );
}

}	// namespace IterationPack
