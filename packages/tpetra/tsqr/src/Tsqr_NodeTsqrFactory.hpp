//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef __TSQR_NodeTsqrFactory_hpp
#define __TSQR_NodeTsqrFactory_hpp

#include <Tsqr_ConfigDefs.hpp>
#include <Kokkos_DefaultNode.hpp>

#ifdef HAVE_KOKKOSTSQR_TBB
#  include <TbbTsqr.hpp>
#endif // HAVE_KOKKOSTSQR_TBB

#include <Tsqr_KokkosNodeTsqr.hpp>
#include <Tsqr_SequentialTsqr.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListExceptions.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include <stdexcept>


namespace TSQR {

  /// \class NodeTsqrFactory
  /// \brief Factory for creating an instance of the right \c NodeTsqr subclass.
  /// \author Mark Hoemmen
  ///
  /// \tparam Node The Kokkos Node type
  /// \tparam Scalar The type of entries in the matrices to factor
  /// \tparam LocalOrdinal The type of local indices in the matrices to factor
  ///
  /// This class maps from a particular Kokkos \c Node type, to the
  /// corresponding \c NodeTsqr subclass.  It lets you construct a
  /// default ParameterList for that \c NodeTsqr subclass, as well as
  /// an instance of the \c NodeTsqr subclass.  It also provides
  /// typedefs for template metaprogramming.
  ///
  /// The "right" \c NodeTsqr subclass is a function of the \c Node
  /// template parameter, and possibly also of the other template
  /// parameters.
  ///
  /// \note If this class does <i>not</i> have a partial
  ///   specialization for your \c Node type, it defaults to use
  ///   SequentialTsqr.  That class does <i>not</i> use threads, and
  ///   only knows how to deal with host data; it cannot handle GPU
  ///   device-resident data.  Thus, it may perform poorly.
  template<class Node, class Scalar, class LocalOrdinal>
  class NodeTsqrFactory {
  public:
    //! The Kokkos Node type.
    typedef Node node_type;
    //! Pointer (RCP) to node_type.
    typedef Teuchos::RCP<node_type> node_ptr;

    //! The NodeTsqr subclass corresponding to the Kokkos Node type.
    typedef SequentialTsqr<LocalOrdinal, Scalar> node_tsqr_type;

    /// \brief Default parameter list for intranode TSQR.
    ///
    /// \note The default implementation returns an empty (not null)
    ///   parameter list.  Each specialization for a specific Node
    ///   type redefines this method to return a parameter list
    ///   appropriate for that Node type's TSQR implementation.
    static Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> params = parameterList ("NodeTsqr");
      // Create a temporary node_tsqr_type instance in order to get
      // default parameters.  The empty input parameter list will get
      // filled in with default values of missing parameters.
      node_tsqr_type nodeTsqr (params);

      return params;
    }

    /// \brief Return a pointer to the intranode TSQR implementation.
    ///
    /// \param node [in/out] Pointer to the Kokkos Node instance.
    ///
    /// \param plist [in/out] Parameter list for configuring the
    ///   NodeTsqr implementation.
    static Teuchos::RCP<node_tsqr_type>
    makeNodeTsqr (const Teuchos::RCP<node_type>& node,
                  const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      (void) node;
      return rcp (new node_tsqr_type (plist));
    }

    /// \brief Prepare the NodeTsqr instance for use by setting its
    ///   Kokkos \c Node instance.
    ///
    /// Some NodeTsqr subclasses can't compute anything until they
    /// have a pointer to a Kokkos Node instance.  Call this method
    /// before invoking any computational methods of the NodeTsqr
    /// subclass instance.
    ///
    /// \pre <tt> ! nodeTsqr.is_null() && ! node.is_null() </tt>
    /// \post <tt> nodeTsqr->ready() </tt>
    static void
    prepareNodeTsqr (const Teuchos::RCP<node_tsqr_type>& nodeTsqr,
                     const Teuchos::RCP<node_type>& node)
    {
      // SequentialTsqr doesn't need the Kokkos Node instance.
      (void) nodeTsqr;
      (void) node;
    }
  };
} // namespace TSQR

#endif // __TSQR_NodeTsqrFactory_hpp
