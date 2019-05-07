// *****************************************************************************
/*!
  \file      src/Inciter/SchemeBase.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Base class to Scheme, a generic interface to discretization proxies
  \details   This file defines the base class to Scheme, a generic interface to
    discretization proxies. This class is intended to be used in conjunction
    with Scheme. See Scheme for usage and more details on how to extend this
    class.
*/
// *****************************************************************************
#ifndef SchemeBase_h
#define SchemeBase_h

#include <tuple>

#include "Variant.h"

#include "Exception.h"
#include "PUPUtil.h"
#include "Inciter/Options/Scheme.h"

#include "NoWarning/diagcg.decl.h"
#include "NoWarning/alecg.decl.h"
#include "NoWarning/distfct.decl.h"
#include "NoWarning/dg.decl.h"
#include "NoWarning/discretization.decl.h"

namespace inciter {

//! Base class for generic forwarding interface to discretization proxies
class SchemeBase {

  private:
    //! Variant type listing all chare proxy types modeling the same concept
    using Proxy =
      std::variant< CProxy_DiagCG, CProxy_DG, CProxy_ALECG >;

  public:
    //! Variant type listing all chare element proxy types (behind operator[])
    using ProxyElem =
      std::variant< CProxy_DiagCG::element_t,
                    CProxy_DG::element_t,
                    CProxy_ALECG::element_t >;

    //! Empty constructor for Charm++
    explicit SchemeBase() {}

    //! Constructor
    //! \param[in] scheme Discretization scheme
    //! \details Based on the input enum we create at least two empty chare
    //!   arrays: (1) discproxy which contains common functionality and data for
    //!   all discretizations, and (2) proxy, which have functionality and data
    //!   specific to a given discretization. Note that proxy is bound (in
    //!   migration behavior and properties) to discproxy.
    //! \note There may be other bound proxy arrays created depending on the
    //!   specific discretization configured by the enum.
    explicit SchemeBase( ctr::SchemeType scheme ) :
      discproxy( CProxy_Discretization::ckNew() )
    {
      bound.bindTo( discproxy );
      if (scheme == ctr::SchemeType::DiagCG) {
        proxy = static_cast< CProxy_DiagCG >( CProxy_DiagCG::ckNew(bound) );
        fctproxy = CProxy_DistFCT::ckNew(bound);
      } else if (scheme == ctr::SchemeType::DG ||
                 scheme == ctr::SchemeType::DGP1 ||
                 scheme == ctr::SchemeType::DGP2 ||
                 scheme == ctr::SchemeType::PDG)
      {
        proxy = static_cast< CProxy_DG >( CProxy_DG::ckNew(bound) );
      } else if (scheme == ctr::SchemeType::ALECG) {
        proxy = static_cast< CProxy_ALECG >( CProxy_ALECG::ckNew(bound) );
      } else Throw( "Unknown discretization scheme" );
    }

    //! Get reference to discretization proxy
    //! \return Discretization Charm++ chare array proxy
    CProxy_Discretization& disc() noexcept { return discproxy; }

    //! Get reference to DistFCT proxy
    //! \return DistFCT Charm++ chare array proxy
    CProxy_DistFCT& fct() noexcept { return fctproxy; }

    //! Get reference to scheme proxy
    //! \return Variant storing Charm++ chare array proxy configured
    const Proxy& getProxy() noexcept { return proxy; }

    //! Query underlying proxy type
    //! \return Zero-based index into the set of types of Proxy
    std::size_t index() const noexcept { return proxy.index(); }

    //! Query underlying proxy element type
    //! \return Zero-based index into the set of types of ProxyElem
    std::size_t index_element() const noexcept {
      return tk::element< ProxyElem >( proxy, 0 ).index();
    }

    //! Charm++ array options accessor for binding external proxies
    //! \return Charm++ array options object reference
    const CkArrayOptions& arrayoptions() { return bound; }

  protected:
    //! Variant storing one proxy to which this class is configured for
    Proxy proxy;
    //! Charm++ proxy to data and code common to all discretizations
    CProxy_Discretization discproxy;
    //! Charm++ proxy to flux-corrected transport (FCT) driver class
    CProxy_DistFCT fctproxy;
    //! Charm++ array options for binding chares
    CkArrayOptions bound;

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      p | proxy;
      p | discproxy;
      p | fctproxy;
      p | bound;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] s SchemeBase object reference
    friend void operator|( PUP::er& p, SchemeBase& s ) { s.pup(p); }
    //@}
};

} // inciter::

#endif // SchemeBase_h
