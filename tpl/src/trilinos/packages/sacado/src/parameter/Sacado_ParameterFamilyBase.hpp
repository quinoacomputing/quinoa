// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_PARAMETERFAMILYBASE_HPP
#define SACADO_PARAMETERFAMILYBASE_HPP

#include <map>
#include <string>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Sacado_mpl_apply.hpp"

namespace Sacado {

  /*!
   * A class to store multiple template instantiations of a single templated
   * parameter.
   */
  template <typename EntryBase, typename EntryType>
  class ParameterFamilyBase  {

  public:

    //! Constructor
    ParameterFamilyBase(const std::string& name,
                        bool supports_ad,
                        bool supports_analytic);

    //! Destructor
    virtual ~ParameterFamilyBase();

    //! Get the name of the family
    std::string getName() const;

    //! Indicates whether parameter supports AD derivatives
    bool supportsAD() const;

    //! Indicates whether parameter supports analytic derivatives
    bool supportsAnalytic() const;

    //! Determine if family has an entry for the given type \c EvalType.
    template <typename EvalType>
    bool hasType() const;

    //! Add a new parameter using custom entry
    /*!
     * Returns true if successful in adding entry to library, false
     * otherwise.  If \c allow_overwrite is true, any existing entry will
     * be overwritten by the supplied entry.
     */
    template <typename EvalType>
    bool
    addEntry(const Teuchos::RCP< typename Sacado::mpl::apply<EntryType,EvalType>::type >& entry,
             const bool allow_overwrite = false);

    //! Gets the entry corresponding to type \em EvalType
    template <typename EvalType>
    Teuchos::RCP< typename Sacado::mpl::apply<EntryType,EvalType>::type >
    getEntry();

    //! Gets the entry corresponding to type \em EvalType
    template <typename EvalType>
    Teuchos::RCP< const typename Sacado::mpl::apply<EntryType,EvalType>::type >
    getEntry() const;

    //! Print the family
    /*!
     * Set print_values = true to print each parameter value for
     * each evaluation type.
     */
    void print(std::ostream& os, bool print_values = false) const;

  protected:

    //! Map of entries for a parameter name
    typedef std::map<std::string, Teuchos::RCP<EntryBase> > EvalMap;

    //! Const iterator for EvalMap
    typedef typename EvalMap::const_iterator const_iterator;

    //! Iterator for EvalMap
    typedef typename EvalMap::iterator iterator;

    //! Returns a string representation of type \em EntryType
    template <class EvalType> std::string getTypeName() const;

  private:

    //! Private to prohibit copying
    ParameterFamilyBase(const ParameterFamilyBase&);

    //! Private to prohibit copying
    ParameterFamilyBase& operator = (const ParameterFamilyBase&);

  protected:

    //! Family of parameter entries
    EvalMap family;

    //! Family name
    const std::string name;

    //! Family supports AD
    bool supports_ad;

    //! Family supports analytic derivatives
    bool supports_analytic;

  };
}

// Include template definitions
#include "Sacado_ParameterFamilyBaseImp.hpp"

#endif
