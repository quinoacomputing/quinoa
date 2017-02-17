/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include <Epetra_CombineMode.h>
#include <ifp_parameters.h>

#ifdef HAVE_TEUCHOS_EXTENDED
#include <Teuchos_StrUtils.hpp>
#endif

namespace Ifpack {

//----------------------------------------------------------------------------
Teuchos::map<string,parameter>& key_map()
{
  static Teuchos::map<string,parameter> ifpack_key_map;
  return( ifpack_key_map );
}

//----------------------------------------------------------------------------
void initialize_string_map()
{
  static bool already_initialized = false;
  if (already_initialized) {
    return;
  }

  Teuchos::map<string,parameter>& ifp_key_map = key_map();

  ifp_key_map["LEVEL_FILL"]    = level_fill;
  ifp_key_map["LEVEL_OVERLAP"] = level_overlap;
  ifp_key_map["ABSOLUTE_THRESHOLD"] = absolute_threshold;
  ifp_key_map["RELATIVE_THRESHOLD"] = relative_threshold;
  ifp_key_map["OVERLAP_MODE"] = overlap_mode;
  ifp_key_map["DROP_TOLERANCE"] = drop_tolerance;
  ifp_key_map["FILL_TOLERANCE"] = fill_tolerance;
  ifp_key_map["RELAX_VALUE"] = relax_value;
  ifp_key_map["USE_RECIPROCAL"] = use_reciprocal;
  ifp_key_map["NUM_STEPS"] = num_steps;

  already_initialized = true;
}

//----------------------------------------------------------------------------
string upper_case(const string& s)
{
#ifdef HAVE_TEUCHOS_EXTENDED
  string upp = Teuchos::StrUtils::allCaps(s);
#else
  string upp(s);
  for(unsigned i=0; i<upp.length(); ++i) {
    upp[i] = toupper(upp[i]);
  }
#endif

  return(upp);
}

//----------------------------------------------------------------------------
void set_parameters(const Teuchos::ParameterList& parameterlist,
                    param_struct& params,
                    bool cerr_warning_if_unused)
{
  initialize_string_map();

  Teuchos::map<string,parameter>& ifp_key_map = key_map();

  Teuchos::ParameterList::ConstIterator
    pl_iter = parameterlist.begin(),
    pl_end  = parameterlist.end();

  for(; pl_iter != pl_end; ++pl_iter) {
    string name = upper_case((*pl_iter).first);

    const Teuchos::ParameterEntry& entry = (*pl_iter).second;
    bool entry_used = false;

    Teuchos::map<string,parameter>::iterator result = ifp_key_map.find(name);
    if (result != ifp_key_map.end()) {
      int dummy_int = -1;
      double dummy_double = -99.9;
      bool dummy_bool = false;
      Epetra_CombineMode dummy_mode = Add;

      parameter offset = (*result).second;

      if (entry.isType<double>()) {
        if (offset < FIRST_INT_PARAM) {
          params.double_params[offset] = entry.getValue(&dummy_double);
          entry_used = true;
        }
      }
      else if (entry.isType<int>()) {
        int int_val = entry.getValue(&dummy_int);
        if (offset >= FIRST_INT_PARAM && offset <= LAST_INT_PARAM) {
          params.int_params[offset-FIRST_INT_PARAM] = int_val;
          entry_used = true;
        }
        else if (offset == use_reciprocal) {
          params.use_reciprocal = int_val;
          entry_used = true;
        }
      }
      else if (entry.isType<bool>()) {
        params.use_reciprocal = entry.getValue(&dummy_bool);
        entry_used = true;
      }
      else if (entry.isType<Epetra_CombineMode>()) {
        params.overlap_mode = entry.getValue(&dummy_mode);
        entry_used = true;
      }
    }

    if (!entry_used && cerr_warning_if_unused) {
      cerr << "Ifpack set_parameters warning: '"<<name<<"' not used."<<endl;
    }
  }
}

} // namespace Ifpack

