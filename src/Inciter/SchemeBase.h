// *****************************************************************************
/*!
  \file      src/Inciter/SchemeBase.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
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

#include "NoWarning/matcg.decl.h"
#include "NoWarning/diagcg.decl.h"
#include "NoWarning/dg.decl.h"
#include "NoWarning/discretization.decl.h"

namespace inciter {

class SchemeBase {

  public:
    //! Empty constructor for Charm++
    explicit SchemeBase() {}

    //! Constructor
    //! \param[in] scheme Discretization scheme
    //! \details Based on the enum we create two empty chare arrays: (1)
    //!    discproxy which contains common functionality and data for all
    //!    discretizations, and (2) proxy, which have functionality and data
    //!    specific to a given discretization. Note that proxy is bound (in
    //!    migration behavior and properties) to discproxy.
    explicit SchemeBase( ctr::SchemeType scheme ) :
      discproxy( CProxy_Discretization::ckNew() )
    {
      CkArrayOptions bound;
      bound.bindTo( discproxy );
      if (scheme == ctr::SchemeType::MatCG) {
        proxy = static_cast< CProxy_MatCG >( CProxy_MatCG::ckNew(bound) );
      } else if (scheme == ctr::SchemeType::DiagCG) {
        proxy = static_cast< CProxy_DiagCG >( CProxy_DiagCG::ckNew(bound) );
      } else if (scheme == ctr::SchemeType::DG) {
        proxy = static_cast< CProxy_DG >( CProxy_DG::ckNew(bound) );
      } else Throw( "Unknown discretization scheme" );
    }

    //! Get reference to discretization proxy
    CProxy_Discretization& get() noexcept { return discproxy; }

    //! Query underlying proxy type
    //! \return Zero-based index into the set of types of Proxy
    int which() const noexcept { return proxy.which(); }

    //! Query underlying proxy element type
    //! \return Zero-based index into the set of types of ProxyElem
    int which_element() const noexcept { return element(0).which(); }

  protected:
    //! Variant type listing all chare proxy types modeling the same concept
    using Proxy = boost::variant< CProxy_MatCG, CProxy_DiagCG, CProxy_DG >;
    //! Variant type listing all chare element proxy types (behind operator[])
    using ProxyElem =
      boost::variant< CProxy_MatCG::element_t, CProxy_DiagCG::element_t,
                      CProxy_DG::element_t >;

    //! Variant storing one proxy to which this class is configured for
    Proxy proxy;
    //! Charm++ proxy to data and code common to all discretizations
    CProxy_Discretization discproxy;

    //! Function dereferencing operator[] on chare proxy returning element proxy
    //! \param[in] x Chare array element index
    //! \return Chare array element proxy as a variant
    //! \details The returning element proxy is a variant, depending on the
    //!   input proxy.
    ProxyElem element( const CkArrayIndex1D& x ) const {
      return boost::apply_visitor( Idx(x), proxy );
    }

    //! Functor to dereference operator[] of chare proxy returning element proxy
    //! \details The returning element proxy is a variant with a type depending
    //!   on the input proxy.
    struct Idx : boost::static_visitor< ProxyElem > {
      Idx( const CkArrayIndex1D& idx ) : x(idx) {}
      template< typename P >
        ProxyElem operator()( const P& p ) const { return p[x]; }
      CkArrayIndex1D x;
    };

    //! Generic base for all call_* classes
    //! \details This class stores the entry method arguments and contains a
    //!   helper function and classes that facilitate unpacking the tuple into a
    //!   variadic list of arguments passed to the called entry method. The Spec
    //!   template argument is the class inheriting from this class (see CRTP at
    //!   https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern),
    //!   specializing functionality of this class. CRTP enables hiding all
    //!   generic functionality here while exposing only the specialized one
    //!   (the actual function call) to the base, which minimizes client code
    //!   (in class Scheme). Unpacking the tuple to a variadic argument list is
    //!   loosely inspired by https://stackoverflow.com/a/16868151.
    template< class Spec, typename... Args >
    struct Call : boost::static_visitor<> {
      //! Constructor storing called member function arguments in tuple
      Call( Args&&... args ) : arg( std::forward_as_tuple(args...) ) {}
      //! Helper class for unpacking tuple into argument list
      template< typename P, typename Tuple, bool Done, int Total, int... N >
      struct invoke_impl {
        static void invoke( P& p, Tuple&& t ) {
          invoke_impl< P, Tuple, Total == 1 + sizeof...(N), Total, N...,
                     sizeof...(N) >::invoke( p, std::forward<Tuple>(t) );
        }
      };
      //! Helper class for unpacking tuple into argument list (end of list)
      template< typename P, typename Tuple, int Total, int... N >
      struct invoke_impl< P, Tuple, true, Total, N... > {
        static void invoke( P& p, Tuple&& t ) {
          Spec::invoke( p, std::get< N >( std::forward<Tuple>(t) )... );
        }
      };
      //! Invoke member function with arguments from tuple
      //! \param[in,out] p Proxy behind which the entry method is called
      //! \param[in] t Optional tuple holding argments to be passed to the call
      template< typename P, typename Tuple = std::tuple<int> >
      static void invoke( P& p, Tuple&& t = {} ) {
        typedef typename std::decay<Tuple>::type ttype;
        invoke_impl< P, Tuple, std::tuple_size<ttype>::value == 0,
          std::tuple_size<ttype>::value >::invoke( p, std::forward<Tuple>(t) );
      }
      //! Function call operator overloading all types used with variant visitor
      //! \param[in,out] p Proxy behind which the entry method is called
      //! \details Classes inheriting from this class should also inherit from
      //!   boost::static_visitor which requires overloading operator(),
      //!   unambiguously accepting any value of types Ts over which classes
      //!   inheriting from this class may be used as functors for static
      //!   visitor with a variant of types Ts.
      template< typename P > void operator()(P& p) const { invoke(p,arg); }
      //! Tuple storing arguments of the entry method to be called
      std::tuple< Args... > arg;
    };

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      auto v = tk::Variant< CProxy_MatCG, CProxy_DiagCG, CProxy_DG >( proxy );
      p | v;
      proxy = v.get();
      p | discproxy;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] s SchemeBase object reference
    friend void operator|( PUP::er& p, SchemeBase& s ) { s.pup(p); }
    //@}
};

} // inciter::

#endif // SchemeBase_h
