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

#include <ostream>
#include <iomanip>
#include <typeinfo>

#include "IterationPack_AlgorithmState.hpp"
#include "Teuchos_Assert.hpp"

namespace {
namespace {
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
} // end namespace
inline void output_spaces(std::ostream& out, int spaces)
{	for(int i = 0; i < spaces; ++i) out << ' '; }
}

namespace IterationPack {

AlgorithmState::iq_id_type AlgorithmState::set_iter_quant(
  const std::string& iq_name, const IQ_ptr& iq)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    iq.get() == NULL, std::invalid_argument
    ,"AlgorithmState::set_iter_quant(...) : The iteration quantity witht the name = \'" << iq_name
    << "\' being inserted has iq.get() == NULL!" );
  iq_id_type new_id = iq_.size();
  std::pair<iq_name_to_id_t::iterator,bool>
    r = iq_name_to_id_.insert(iq_name_to_id_t::value_type(iq_name,new_id));
  TEUCHOS_TEST_FOR_EXCEPTION(
    !r.second // an insert did not take place, key = iq_name already existed.
    ,AlreadyExists
    ,"AlgorithmState::set_iter_quant(...) : An iteration quantity with the name \""
    << iq_name << "\" already exists with the iq_id = " << (*r.first).second );
  iq_.push_back(iq);
  return new_id;
}

void AlgorithmState::erase_iter_quant(const std::string& iq_name) {
  iq_name_to_id_t::iterator itr = find_and_assert(iq_name);
  const iq_id_type iq_id = (*itr).second;
  iq_[iq_id] = Teuchos::null;  // set the pointer to null
  iq_name_to_id_.erase( itr );
}

IterQuantity& AlgorithmState::iter_quant(iq_id_type iq_id) {
  bool exists = true;
  try {
    IQ_ptr &_iq = iq_.at(iq_id);
    if( _iq.get() )
      return *_iq;
    else
      exists = false;
  }
  catch(const std::out_of_range& excpt) {	// Thrown by MS VC++ 6.0
    exists = false;
  }
  catch(const std::range_error& excpt) {	// Thrown by libstdc++ v3 in g++ 2.95.2
    exists = false;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    !exists, DoesNotExist
    ,"AlgorithmState::iter_quant(iq_id) : Error, the iteration quantity iq_id = "
    << iq_id << " does not exist.  "
    << ( iq_id < iq_.size()
       ? "This iteration quantity was set and then erased."
       : "This iteration quantity was never set by the client." ) );
  return *iq_.at(0);	// Will never be executed.
}

const IterQuantity& AlgorithmState::iter_quant(iq_id_type iq_id) const {
  return const_cast<AlgorithmState*>(this)->iter_quant(iq_id);
}

void AlgorithmState::next_iteration(bool incr_k) {
  if(incr_k) this->incr_k();
  for(iq_t::iterator itr = iq_.begin(); itr != iq_.end(); ++itr)
    if(itr->get()) (*itr)->next_iteration();
}

void AlgorithmState::dump_iter_quant(std::ostream& out) const {
  using std::setw;
  using std::endl;

  // Find the maximum length of an iteration quantity name
  int name_w_max = 0;
  {for(	iq_name_to_id_t::const_iterator itr = iq_name_to_id_.begin();
      itr !=  iq_name_to_id_.end(); ++itr )
  {
    name_w_max = my_max( (int)name_w_max, (int)(*itr).first.length() );
  }}	

  const int name_w = name_w_max + 4, id_w = 6;
  const char gap[] = "    ";

  out		<< "\n\n"
      << std::left    << setw(name_w) << "iq_name"
      << std::right   << setw(id_w)   << "iq_id"
      << gap          << std::left    << "concrete type of iq / concrete type of object\n";

  out		<< std::left    << setw(name_w) << "-----"
      << std::right   << setw(id_w)   << "------"
      << gap          << std::left    << "---------------------------------------------\n";

  {for(	iq_name_to_id_t::const_iterator itr = iq_name_to_id_.begin();
      itr !=  iq_name_to_id_.end(); ++itr )
  {
    out << std::left    << setw(name_w) << (*itr).first
      << std::right   << setw(id_w)   << (*itr).second
      << gap          << std::left    << typeName(*iq_[(*itr).second]) << endl
      << std::left    << setw(name_w) << ""
      << std::right   << setw(id_w)   << ""
      << gap          << std::left;
    iq_[(*itr).second]->print_concrete_type(out);
    out << endl;
  }}
}

// private

AlgorithmState::iq_name_to_id_t::iterator AlgorithmState::find_and_assert(
  const std::string& iq_name)
{
  iq_name_to_id_t::iterator itr = iq_name_to_id_.find(iq_name);
  if(itr == iq_name_to_id_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, DoesNotExist
      ,"AlgorithmState::find_and_assert(iq_name) : The iteration "
      "quantity with the name \"" << iq_name << "\" does not exist" );
  return itr;
}

AlgorithmState::iq_name_to_id_t::const_iterator AlgorithmState::find_and_assert(
  const std::string& iq_name) const
{
  iq_name_to_id_t::const_iterator itr = iq_name_to_id_.find(iq_name);
  if(itr == iq_name_to_id_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, DoesNotExist
      ,"AlgorithmState::find_and_assert(iq_name) : The iteration "
      "quantity with the name \"" << iq_name << "\" does not exist" );
  return itr;
}

}	// end namespace IterationPack
