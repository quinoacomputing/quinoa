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
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>

#include "ProfileHackPack_profile_hack.hpp"

namespace {

//
class TimingEntry {
public:
  TimingEntry() : num_calls(0), total_time_secs(0) {}
  size_t         num_calls;
  double         total_time_secs;   // in seconds
}; // end struct TimingEntry

//
typedef std::map<std::string,TimingEntry>  func_timing_map_t;

//
func_timing_map_t   func_timing_map;

//
class SortByTimingDescending {
public:
  bool operator()( const func_timing_map_t::value_type& x
    , const func_timing_map_t::value_type& y ) const
  {
    return x.second.total_time_secs > y.second.total_time_secs;
  }
};	// end class AbsMultVal

//
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }

} // end namespace 

void ProfileHackPack::set_time( const char func_name[], double time_secs )
{
  TimingEntry &entry = func_timing_map[func_name];
  entry.num_calls++;
  entry.total_time_secs += time_secs;
}

void ProfileHackPack::print_timings( std::ostream& out )
{
  using std::setw;
  using std::right;
  using std::left;

  // Sort the entries by the function times in descending order
  typedef std::vector<std::pair<std::string,TimingEntry> > list_sorted_t;
  list_sorted_t  list_sorted(func_timing_map.size());
  {
    func_timing_map_t::const_iterator itr_from = func_timing_map.begin();
    list_sorted_t::iterator           itr_to   = list_sorted.begin();
    for( ; itr_from != func_timing_map.end(); ++itr_from, ++itr_to ) {
      itr_to->first = itr_from->first;
      itr_to->second = itr_from->second;
    }
  }
  std::copy( func_timing_map.begin(), func_timing_map.end(), list_sorted.begin() );
  std::sort( list_sorted.begin(), list_sorted.end(), SortByTimingDescending() );
  // Get the maximum function name size
  int max_func_name_len = 25;
  {for( list_sorted_t::const_iterator itr = list_sorted.begin(); itr != list_sorted.end(); ++itr )
    max_func_name_len = my_max( int(max_func_name_len), int(itr->first.size()) );}
  // Print out the function times
  const int
    name_w = max_func_name_len+2,
    dbl_w  = 22,
    int_w  = 10;
  const char
    name_ul[] = "-------------------------",
    dbl_ul[]  = "--------------------",
    int_ul[]  = "--------";

  out << "\nPoor man\'s profile times:\n\n";
  out << left  << setw(name_w) << "function name"
    << right << setw(dbl_w)  << "self+childern(sec)"
    << right << setw(int_w)  << "# calls"
    << right << setw(dbl_w)  << "av cpu/call(sec)"
    << std::endl
    << left  << setw(name_w) << name_ul
    << right << setw(dbl_w)  << dbl_ul
    << right << setw(int_w)  << int_ul
    << right << setw(dbl_w)  << dbl_ul
    << std::endl;
  {for( list_sorted_t::const_iterator itr = list_sorted.begin(); itr != list_sorted.end(); ++itr ) {
  out << left  << setw(name_w) << itr->first
    << right << setw(dbl_w)  << itr->second.total_time_secs
    << right << setw(int_w)  << itr->second.num_calls
    << right << setw(dbl_w)  << (itr->second.total_time_secs / itr->second.num_calls)
    << std::endl;
  }}
}
