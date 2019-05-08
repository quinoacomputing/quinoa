// *****************************************************************************
/*!
  \file      src/Base/Variant.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Helpers for operator operator[] using boost::variant
  \details   Helpers for applying operator[] using boost::variant.
*/
// *****************************************************************************
#ifndef Variant_h
#define Variant_h

#include "NoWarning/variant.hpp"
#include "NoWarning/charm++.hpp"

namespace tk {

//! Functor to dereference operator[] of chare proxy inside a variant
//! \details Since the chare array proxy is behind a variant, the returning
//!   element proxy from operator() is also a variant, defined by ProxyElem with
//!   a type depending on the input proxy, given by P, i.e., overloaded for all
//!   proxy types the variant supports.
template< class ProxyElem >
struct Idx : boost::static_visitor< ProxyElem > {
  explicit Idx( const CkArrayIndex1D& idx ) : x(idx) {}
  template< typename P >
    ProxyElem operator()( const P& p ) const { return p[x]; }
  CkArrayIndex1D x;
};

//! Function dereferencing operator[] of chare proxy inside variant
//! \param[in] proxy Chare array proxy inside a variant to dereference
//! \param[in] x Chare array element index
//! \return Chare array element proxy as a variant, defined by ProxyElem
//! \details The returning element proxy is a variant, depending on the input
//!   proxy.
//! \see inciter::Scheme, inciter::SchemeBase, or, e.g., DistFCT::apply() for
//!   client code.
template< class ProxyElem, class Proxy >
ProxyElem element( const Proxy& proxy, const CkArrayIndex1D& x ) {
  return boost::apply_visitor( Idx<ProxyElem>(x), proxy );
}

} // tk::

#endif // Variant_h
