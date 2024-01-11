// *****************************************************************************
/*!
  \file      src/Inciter/Scheme.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Polymorphic glue for calling Charm++ entry methods to base class
    Discretization, its children implementing specific discretization schemes,
    and helper classes
  \details
    The purpose of this class is to hide, behind a single type, different
    Charm++  proxy types that model a single concept, i.e., define some common
    functions as Charm++ entry methods that can be used in either a broadcast
    and/or in a way of addressing a single array element. As a result, member
    functions can be invoked by client code without knowing the underlying type
    or any specifics to the underlying differences of the classes that model the
    same concept, i.e., expose the same member functions. The idea is very
    similar to inheritance and runtime polymorphism with base classes and
    virtual functions: some member functions and data are common to all types
    modeled (and thus are not repeated and/or copied), while some are specific.
    A difference is that the "base" and "child" classes are Charm++ proxies.
    Note that while Charm++ does support inheritance and runtime polymorphism
    with chare arrays, we still prefer the implementation below because it uses
    entirely value semantics (inside and in client code) and thus it keeps the
    complexity of the dispatch behind this class and does not expose it to
    client code.

    The advantages of this class over traditional runtime polymorphism are (1)
    value semantics (both internally and to client code), (2) not templated,
    and (3) PUPable, i.e., an instance of Scheme can be sent across the network
    using Charm++'s pup framework. Also, since the class only holds a couple of
    chare proxies, it is lightweight.

    Example usage from client code:

    \code{.cpp}
      // Instantiate a Scheme object
      Scheme s( ctr::SchemeType::DG );  // see Control/Inciter/Options/Scheme.h

      // Issue broadcast to child scheme entry method
      s.bcast< Scheme::setup >(...);

      // Issue broadcast to base (Discretization) entry method
      s.disc().totalvol();
    \endcode

    Organization, implementation details, end extension of the class:

    Scheme contains (at least) two Charm++ proxies: discproxy and proxy. The
    former contains data and functionality common to all discretizations, and
    this can be considered as an equivalent to a base class in the OOP sense.
    The latter, proxy, contains data and functionality specific to a particular
    discretization. When instantiated, Scheme is configured for a single
    specific discretization which must be selected from the list of types in
    SchemeBase::Proxy.

    The underlying type of proxy is a variant, which allows storing exactly one
    object. A variant is a type-safe union. An instance of a variant at any
    given time either holds a value of one of its alternative types. Read more
    on std::variant on how they work.

    Adding a new child scheme is done by
    (1) Adding a new type of Charm++ chare array proxy to Scheme::Proxy,
    (2) Adding a new type of Charm++ chare array element proxy to
        Scheme::ProxyElem, and
    (3) Adding a new branch to the if test in Scheme's constructor.

  \see A talk on "Concept-based runtime polymorphism with Charm++ chare arrays
    using value semantics given by J. Bakosi at the 16th Annual Workshop on
    Charm++ and its Applications, April 2018, discussing an earlier, more
    verbose implementation of the idea, using C++11.
*/
// *****************************************************************************
#ifndef Scheme_h
#define Scheme_h

#include "Exception.hpp"
#include "PUPUtil.hpp"
#include "Inciter/Options/Scheme.hpp"

#include "NoWarning/discretization.decl.h"
#include "NoWarning/diagcg.decl.h"
#include "NoWarning/alecg.decl.h"
#include "NoWarning/oversetfe.decl.h"
#include "NoWarning/distfct.decl.h"
#include "NoWarning/dg.decl.h"
#include "NoWarning/fv.decl.h"
#include "NoWarning/ale.decl.h"
#include "NoWarning/conjugategradients.decl.h"
#include "NoWarning/ghosts.decl.h"

namespace inciter {

//! Base class for generic forwarding interface to discretization proxies
class Scheme {

  private:
    //! Variant type listing all chare proxy types modeling the same concept
    using Proxy = std::variant< CProxy_DiagCG
                              , CProxy_DG
                              , CProxy_ALECG
                              , CProxy_OversetFE
                              , CProxy_FV >;

  public:
    //! Variant type listing all chare element proxy types
    using ProxyElem = std::variant< CProxy_DiagCG::element_t
                                  , CProxy_DG::element_t
                                  , CProxy_ALECG::element_t
                                  , CProxy_OversetFE::element_t
                                  , CProxy_FV::element_t >;

    //! Empty constructor for Charm++
    explicit Scheme() {}

    //! Constructor
    //! \param[in] scheme Discretization scheme
    //! \param[in] ale True if enable ALE
    //! \param[in] linearsolver True if enable a linear solver
    //! \details Based on the input enum we create at least two empty chare
    //!   arrays: (1) discproxy which contains common functionality and data for
    //!   all discretizations, and (2) proxy, which have functionality and data
    //!   specific to a given discretization. Note that proxy is bound (in
    //!   migration behavior and properties) to discproxy.
    //! \note There may be other bound proxy arrays created depending on the
    //!   specific discretization configured by the enum.
    explicit Scheme( ctr::SchemeType scheme,
                     bool ale = false,
                     bool linearsolver = false,
                     tk::Centering centering = tk::Centering::NODE ) :
      discproxy( CProxy_Discretization::ckNew() )
    {
      bound.bindTo( discproxy );
      if (scheme == ctr::SchemeType::DiagCG) {
        proxy = static_cast< CProxy_DiagCG >( CProxy_DiagCG::ckNew(bound) );
        fctproxy = CProxy_DistFCT::ckNew(bound);
      } else if (scheme == ctr::SchemeType::DG ||
                 scheme == ctr::SchemeType::P0P1 ||
                 scheme == ctr::SchemeType::DGP1 ||
                 scheme == ctr::SchemeType::DGP2 ||
                 scheme == ctr::SchemeType::PDG)
      {
        proxy = static_cast< CProxy_DG >( CProxy_DG::ckNew(bound) );
      } else if (scheme == ctr::SchemeType::ALECG) {
        proxy = static_cast< CProxy_ALECG >( CProxy_ALECG::ckNew(bound) );
      } else if (scheme == ctr::SchemeType::OversetFE) {
        proxy = static_cast< CProxy_OversetFE >( CProxy_OversetFE::ckNew(bound) );
      } else if (scheme == ctr::SchemeType::FV) {
        proxy = static_cast< CProxy_FV >( CProxy_FV::ckNew(bound) );
      } else Throw( "Unknown discretization scheme" );
      if (ale) aleproxy = CProxy_ALE::ckNew(bound);
      if (linearsolver)
        conjugategradientsproxy = tk::CProxy_ConjugateGradients::ckNew(bound);
      if (centering == tk::Centering::ELEM)
        ghostsproxy = CProxy_Ghosts::ckNew(bound);
    }

    //! Entry method tags for specific Scheme classes to use with bcast()
    struct setup {};
    struct box {};
    struct transferSol {};
    struct advance {};
    struct resized {};
    struct resizeComm {};
    struct refine {};
    struct lhs {};
    struct nodeNeighSetup {};
    struct diag {};
    struct evalLB {};
    struct doneInserting {};
    //! Issue broadcast to Scheme entry method
    //! \tparam Fn Function tag identifying the entry method to call
    //! \tparam Args Types of arguments to pass to entry method
    //! \param[in] args Arguments to member function entry method to be called
    //! \details This function issues a broadcast to a member function entry
    //!   method of the Scheme chare array (the child of Discretization) and is
    //!   thus equivalent to proxy.Fn(...).
    template< typename Fn, typename... Args >
    void bcast( Args&&... args ) {
      std::visit( [&]( auto& p ){
          if constexpr( std::is_same_v< Fn, setup > )
            p.setup( std::forward< Args >( args )... );
          if constexpr( std::is_same_v< Fn, box > )
            p.box( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, transferSol > )
            p.transferSol( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, advance > )
            p.advance( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, resized > )
            p.resized( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, resizeComm > )
            p.resizeComm( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, refine > )
            p.refine( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, lhs > )
            p.lhs( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, nodeNeighSetup > )
            p.nodeNeighSetup( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, diag > )
            p.diag( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, evalLB > )
            p.evalLB( std::forward< Args >( args )... );
          else if constexpr( std::is_same_v< Fn, doneInserting > )
            p.doneInserting( std::forward< Args >( args )... );
        }, proxy );
    }

    //! Function tags for specific Scheme classes to use with ckLocal()
    struct resizePostAMR {};
    struct extractFieldOutput {};
    struct solution {};
    //! Call Scheme function via Charm++ chare array element's ckLocal()
    //! \tparam Fn Function tag identifying the function to call
    //! \tparam Args Types of arguments to pass to function
    //! \param[in] x Chare array element index
    //! \param[in] args Arguments to member function function to be called
    //! \details This function calls a member function via Charm++'s ckLocal()
    //!    behind the element proxy configured, indexed by the array index x.
    //!    Since the call is behind ckLocal(), the member function does not have
    //!    to be a Charm++ entry method.
    template< typename Fn, typename... Args >
    auto ckLocal( const CkArrayIndex1D& x, Args&&... args ) const {
      auto e = element( x );
      return std::visit( [&]( auto& p ){
          if constexpr( std::is_same_v< Fn, resizePostAMR > )
            return p.ckLocal()->resizePostAMR( std::forward<Args>(args)... );
          else if constexpr( std::is_same_v< Fn, extractFieldOutput > )
            return p.ckLocal()->extractFieldOutput(std::forward<Args>(args)...);
          else if constexpr( std::is_same_v< Fn, solution > )
            return p.ckLocal()->solution( std::forward<Args>(args)... );
        }, e );
    }

    //! Function to call the insert entry method of an element proxy
    //! \param[in] x Chare array element index
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the insert member function of a chare array
    //!   element proxy and thus equivalent to proxy[x].insert(...), using the
    //!   last argument as default.
    template< typename... Args >
    void insert( const CkArrayIndex1D& x, Args&&... args ) {
      auto e = element( x );
      std::visit( [&]( auto& p ){ p.insert(std::forward<Args>(args)...); }, e );
    }

    //! Get reference to discretization proxy
    //! \return Discretization Charm++ chare array proxy
    CProxy_Discretization& disc() noexcept { return discproxy; }

    //! Get reference to DistFCT proxy
    //! \return DistFCT Charm++ chare array proxy
    CProxy_DistFCT& fct() noexcept { return fctproxy; }

    //! Get reference to ALE proxy
    //! \return ALE Charm++ chare array proxy
    CProxy_ALE& ale() noexcept { return aleproxy; }

    //! Get reference to ConjugateGradients proxy
    //! \return ConjugateGradients Charm++ chare array proxy
    tk::CProxy_ConjugateGradients& conjugategradients() noexcept
    { return conjugategradientsproxy; }

    //! Get reference to Ghosts proxy
    //! \return Ghosts Charm++ chare array proxy
    CProxy_Ghosts& ghosts() noexcept { return ghostsproxy; }

    //! Get reference to scheme proxy
    //! Get reference to scheme proxy
    //! \return Variant storing Charm++ chare array proxy configured
    const Proxy& getProxy() noexcept { return proxy; }

    //! Query underlying proxy type
    //! \return Zero-based index into the set of types of Proxy
    std::size_t index() const noexcept { return proxy.index(); }

    //! Query underlying proxy element type
    //! \return Zero-based index that can be used, e.g., indexing into the set
    //!   of types of ProxyElem
    std::size_t index_element() const noexcept { return element(0).index(); }

    //! Charm++ array options accessor for binding external proxies
    //! \return Charm++ array options object reference
    const CkArrayOptions& arrayoptions() { return bound; }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      p | proxy;
      p | discproxy;
      p | fctproxy;
      p | aleproxy;
      p | conjugategradientsproxy;
      p | ghostsproxy;
      p | bound;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] s Scheme object reference
    friend void operator|( PUP::er& p, Scheme& s ) { s.pup(p); }
    //@}

  private:
    //! Variant storing one proxy to which this class is configured for
    Proxy proxy;
    //! Charm++ proxy to data and code common to all discretizations
    CProxy_Discretization discproxy;
    //! Charm++ proxy to flux-corrected transport (FCT) driver class
    CProxy_DistFCT fctproxy;
    //! Charm++ proxy to ALE class
    CProxy_ALE aleproxy;
    //! Charm++ proxy to conjugate gradients linear solver class
    tk::CProxy_ConjugateGradients conjugategradientsproxy;
    //! Charm++ proxy to Ghosts class
    CProxy_Ghosts ghostsproxy;
    //! Charm++ array options for binding chares
    CkArrayOptions bound;

    //! Function dereferencing operator[] of chare proxy inside variant
    //! \param[in] x Chare array element index
    //! \return Chare array element proxy as a variant, defined by ProxyElem
    //! \details The returning element proxy is a variant, depending on the
    //!   input proxy.
    ProxyElem element( const CkArrayIndex1D& x ) const {
      return std::visit( [&]( const auto& p ){
               return static_cast< ProxyElem >( p[x] ); }, proxy );
    }
};

} // inciter::

#endif // Scheme_h
