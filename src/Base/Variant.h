// *****************************************************************************
/*!
  \file      src/Base/Variant.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Variant is a wrapper around boost::variant to make it PUPable
  \details   Variant is a wrapper around boost::variant to make it PUPable.

    Since boost::variant (as well as std::variant) when default-constructed is
    initialized to hold a value of the first alternative of its type list,
    calling PUP that works based on a boost::visitor with a templated operator()
    always incorrectly triggers the overload for the first type. Thus when
    PUPing a variant not only its value but its type must also be sent during
    migration. The class below achieves this by reading out not only the value
    but also its zero-based index of the type alternative that is currently held
    by the variant passed to its initializer constructor. The index and the
    value are stored as an integer and a value in a tuple that can hold a value
    of any of the alternative types the variant can. The index and tuple are
    then PUPed and when unpacking, as an additional step, the variant is
    reconstructed using the index and the value in the tuple. This latter is
    done by invoking a runtime loop across all alternative types using
    boost::mpl::for_each. The final variant can then be accessed by the get()
    public member function.

  \see UnitTest/tests/Base/TestPUPUtil.h or Inciter::SchemeBase.h for
    tk::Variant in action.
*/
// *****************************************************************************
#ifndef Variant_h
#define Variant_h

#include <tuple>

#include "NoWarning/variant.h"
#include "NoWarning/for_each.h"
#include "NoWarning/vector.h"
#include "NoWarning/pup.h"
#include "NoWarning/charm++.h"

#include "Get_Tuple.h"

namespace tk {

template< typename... Types >
class Variant {

  public:
    //! Default constructor for Charm++ (e.g., PUP's receive side)
    Variant() {}

    //! Initializer constructor
    //! \param[in] v boost::variant used as input/output
    //! \details This constructor queries the index at which the variant holds a
    //!   value and sets that value in tuple mathing the type.
    Variant( boost::variant< Types... >& v ) : idx( v.which() ), variant(v) {
      boost::apply_visitor( getval(this), v );
    }

    //! Access boost::variant
    //! \return boost::variant held
    boost::variant< Types... > get() { return variant; }

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \details All of Charm++ PUP modes of sizing, packing, and unpacking PUPs
    //!   the index and the tuple. Additionally, when we unpack the variant is
    //!   also set using the index and the value from tuple at position idx.
    void pup( PUP::er &p ) {
      p | idx;
      p | tuple;
      if (p.isUnpacking())
        boost::mpl::for_each< boost::mpl::vector<Types...> >( setval(this) );
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] v Variant object reference
    friend void operator|( PUP::er& p, Variant< Types... >&& v ) { v.pup(p); }
    //@}

  private:
    //! Visitor setting a value of tuple that matches the type of the variant
    struct getval : boost::static_visitor<> {
      Variant* const host;
      getval( Variant* const h ) : host(h) {}
      template< typename P > void operator()( const P& p ) const {
        tk::get< P >( host->tuple ) = p;
      }
    };

    //! Functor setting the variant based on idx
    struct setval {
      Variant* const host;
      int cnt;
      setval( Variant* const h ) : host(h), cnt(0) {}
      template< typename U > void operator()( U ) {
        if (host->idx == cnt++) host->variant = tk::get< U >( host->tuple );
      }
    };

    //! Index at which the variant holds a value
    int idx;
    //! Tuple that can hold any value of the variant
    std::tuple< Types... > tuple;
    //! Input/output variant
    boost::variant< Types... > variant;
};

//! Functor to dereference operator[] of chare proxy inside a variant
//! \details Since the chare array proxy is behind a variant, the returning
//!   element proxy from operator() is also a variant, defined by ProxyElem with
//!   a type depending on the input proxy, given by P, i.e., overloaded for all
//!   proxy types the variant supports.
template< class ProxyElem >
struct Idx : boost::static_visitor< ProxyElem > {
  Idx( const CkArrayIndex1D& idx ) : x(idx) {}
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
