// *****************************************************************************
/*!
  \file      src/Inciter/SchemeBase.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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
      boost::variant< CProxy_DiagCG, CProxy_DG, CProxy_ALECG >;

  public:
    //! Variant type listing all chare element proxy types (behind operator[])
    using ProxyElem =
      boost::variant< CProxy_DiagCG::element_t,
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
                 scheme == ctr::SchemeType::DGP2)
      {
        proxy = static_cast< CProxy_DG >( CProxy_DG::ckNew(bound) );
      } else if (scheme == ctr::SchemeType::ALECG) {
        proxy = static_cast< CProxy_ALECG >( CProxy_ALECG::ckNew(bound) );
      } else Throw( "Unknown discretization scheme" );
    }

    //! Get reference to discretization proxy
    //! \return Discretization Charm++ chare array proxy
    CProxy_Discretization& get() noexcept { return discproxy; }

    //! Get reference to scheme proxy
    //! \return Variant storing Charm++ chare array proxy configured
    const Proxy& getProxy() noexcept { return proxy; }

    //! Query underlying proxy type
    //! \return Zero-based index into the set of types of Proxy
    int which() const noexcept { return proxy.which(); }

    //! Query underlying proxy element type
    //! \return Zero-based index into the set of types of ProxyElem
    int which_element() const noexcept {
      return tk::element< ProxyElem >( proxy, 0 ).which();
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
      explicit Call( Args&&... args ) : arg( std::forward_as_tuple(args...) ) {}
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
