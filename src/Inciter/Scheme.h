// *****************************************************************************
/*!
  \file      src/Inciter/Scheme.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Generic forwarding interface to discretization proxies
  \details   This file defines a generic interface to discretization proxies.

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
    complexity of the dispatch behind this class and does not expose it client
    code.

    The advantages of this class over traditional runtime polymorphism are (1)
    value semantics (both internally and to client code), (2) not templated,
    and (3) PUPable, i.e., an instance of Scheme can be sent across the network
    using Charm++'s pup framework. Also, since the class only holds a couple of
    chare proxies, it is extremely lightweight.

    Example usage from client code:

    \code{.cpp}
      // Instantiate a Scheme object
      Scheme s( ctr::SchemeType::DG );  // see Control/Inciter/Options/Scheme.h

      // Call a member function entry method in broadcast fashion
      s.coord< tag::bcast >( ... );     // equivalent to proxy.coord( ... );

      // Call a member function entry method in addressing a single array
      // element
      s.coord< tag::elem >( 0, ... );   // equivalent to proxy[0].coord( ... );

      // Broadcast to a member function with optinoal CkEntryOptions object
      CkEntryOptions opt;
      s.coord< tag::bcast >( ..., opt );     // proxy.coord( ..., opt );

      // Address array element with optinoal CkEntryOptions object
      s.coord< tag::elem >( 0, ..., opt );   // proxy[0].coord( ..., opt );
    \endcode

    Organization, implementation details, end extension of the class:

    Scheme, via inheriting from SchemeBase, contains two Charm++ proxies:
    discproxy and proxy. The former contains data and functionality common to
    all discretizations, and this can be considered as an equivalent to a base
    class in the OOP sense. The latter, proxy, contains data and functionality
    specific to a particular discretization. When instantiated, Scheme is
    configured for a single specific discretization which must be selected from
    the list of types in SchemeBase::Proxy.

    The underlying type of proxy is a variant, which allows storing exactly one
    object. A variant is a type-safe union. An instance of a variant at any
    given time either holds a value of one of its alternative types. Read more
    on std::variant or boost::variant on how they work.

    All new member functions that comprise of the concept of the underlying
    proxies, i.e., the interface, must be defined in Scheme. Whereas common
    data, functionality, as well as the list of the proxy types that can be
    configured are defined in SchemeBase. Adding a new forwarding function
    either as a broadcast or addressing a particular chare array element can be
    done by simply copying an existing (similar) one and modifying what
    underlying function (entry method) it calls. The ones that forward to
    discproxy are grouped first, while the ones that forward to the specific
    proxy are listed as second. Using SFINAE, multiple overloads are (and can
    be) defined for a single function, depending on (1) whether it is a
    broadcast or addressing an array element, (2) whether it takes an optional
    (default) last argument, usually used for passing a CkEntryOptions object.
    You can see the Charm++-generated .decl.h file to see what (if any) default
    arguments a particular entry method may take.

    Currently, forwarding functions are defined for two types entry method
    calls: broadcasts, i.e., proxy.fn(), and addressing a particular element,
    i.e., proxy[x].fn(). Another third might be useful to add in the future and
    that is addressing an entry method behind a section proxy. As we have not
    used section proxies so far, this is not yet implemented, but if necessary,
    it should be relatively straightforward to do.

    Extending this class to other discretization schemes is done entirely in
    SchemeBase. Adding a new discretization scheme amounts to, at the minimum:
    (1) Adding a new type of Charm++ chare array proxy to SchemeBase::Proxy,
    (2) Adding a new type of Charm++ chare array element proxy to
        SchemeBase::ProxyElem, and
    (3) Adding a new branch to the if test in SchemeBase's constructor.

    Implementation details: All forwarding calls are implemented taking a
    variadic parameter pack, which can take any number of arguments (including
    zero) and use perfect forwarding to pass the arguments to the entry method
    calls. This way the code remains generic, easy to modify, and the compiler
    automatically adjusts the generated forwarding calls if the types and/or
    number of the arguments to any of the entry methods change. One exception to
    this is those forwarding calls that deal with default arguments, allowing
    for passing CkEntryOptions objects. There the number of arguments are
    hard-coded in the SFINAE construct, but should also be straightforward to
    modify if necessary.

    The functors named as call_* are used to dispatch (at compile time) entry
    method calls behind proxy, whose type is different depending on what
    specific discretization type is configured in the constructor. All common
    functionality in the call_* functors are lifted over to SchemeBase::Call.
    This helps keeping the function-call-specific code in Scheme minimal and
    reuses the generic part in SchemeBase.

    Note that another way of doing the dispatch, that is now done using the
    call_* functors, could have been implemented using a (compile-, or runtime)
    associative container storing std::function objects which would store
    pre-bound function arguments. That would work, but there are three problems
    with such an approach: (1) std::function is not obvious how to pup, i.e.,
    send across the network, (2) std::bind cannot currently be used to bind a
    variadic number arguments and thus the bind calls would not be very generic,
    and (3) a runtime associative container would take additional state.
*/
// *****************************************************************************
#ifndef Scheme_h
#define Scheme_h

#include "Tags.h"
#include "SchemeBase.h"

namespace inciter {

class Scheme : public SchemeBase {

  public:
    // Inherit base constructors
    using SchemeBase::SchemeBase;

    // Calls to discproxy, common to all discretizations

    //////  discproxy.coord(...)
    //! \brief Function to call the coord() entry method of an array discproxy
    //!   (broadcast)
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the coord member function of a chare array
    //!   discproxy and thus equivalent to discproxy.coord(...).
    template< class Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::bcast >::value, int >::type = 0 >
    void coord( Args&&... args ) {
      discproxy.coord( std::forward<Args>(args)... );
    }
    //////  discproxy[x].coord(...)
    //! Function to call the coord() entry method of an element discproxy (p2p)
    //! \param[in] x Chare array element index
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the coord member function of a chare array
    //!   element discproxy and thus equivalent to discproxy[x].coord(...).
    template< typename Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::elem >::value, int >::type = 0 >
    void coord( const CkArrayIndex1D& x, Args&&... args ) {
      discproxy[x].coord( std::forward<Args>(args)... );
    }

    //////  discproxy.totalvol(...)
    //! \brief Function to call the totalvol() entry method of an array
    //!   discproxy (broadcast)
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the totalvol member function of a chare
    //!   array discproxy and thus equivalent to discproxy.totalvol(...).
    template< class Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::bcast >::value, int >::type = 0 >
    void totalvol( Args&&... args ) {
      discproxy.totalvol( std::forward<Args>(args)... );
    }
    //////  discproxy[x].totalvol(...)
    //! \brief Function to call the totalvol() entry method of an element
    //!    discproxy (p2p)
    //! \param[in] x Chare array element index
    //! \param[in] args Arguments to member function (entry method) to be
    //!    called
    //! \details This function calls the totalvol member function of a chare
    //!   array element discproxy and thus equivalent to
    //!   discproxy[x].totalvol(...).
    template< typename Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::elem >::value, int >::type = 0 >
    void totalvol( const CkArrayIndex1D& x, Args&&... args ) {
      discproxy[x].totalvol( std::forward<Args>(args)... );
    }

    //////  discproxy.stat(...)
    //! \brief Function to call the stat() entry method of an array discproxy
    //!   (broadcast)
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the stat member function of a chare array
    //!   discproxy and thus equivalent to discproxy.stat(...).
    template< class Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::bcast >::value, int >::type = 0 >
    void stat( Args&&... args ) {
      discproxy.stat( std::forward<Args>(args)... );
    }
    //////  discproxy[x].stat(...)
    //! Function to call the stat() entry method of an element discproxy (p2p)
    //! \param[in] x Chare array element index
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the stat member function of a chare array
    //!   element discproxy and thus equivalent to discproxy[x].stat(...).
    template< typename Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::elem >::value, int >::type = 0 >
    void stat( const CkArrayIndex1D& x, Args&&... args ) {
      discproxy[x].stat( std::forward<Args>(args)... );
    }

    //////  discproxy[x].insert(...)
    //! Function to call the insert entry method of an element discproxy (p2p)
    //! \param[in] x Chare array element index
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the insert member function of a chare array
    //!   element discproxy and thus equivalent to discproxy[x].insert(...),
    //!   using the last argument as default.
    template< typename... Args >
    void discInsert( const CkArrayIndex1D& x, Args&&... args ) {
      discproxy[x].insert( fctproxy, std::forward<Args>(args)... );
    }

    //////  discproxy.doneInserting(...)
    //! \brief Function to call the doneInserting entry method of an array
    //!   discproxy (broadcast)
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the doneInserting member function of a
    //!   chare array discproxy and thus equivalent to
    //!   discproxy.doneInserting(...).
    template< class Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::bcast >::value, int >::type = 0 >
    void doneDiscInserting( Args&&... args ) {
      discproxy.doneInserting( std::forward<Args>(args)... );
    }

    //////  fctcproxy.doneInserting(...)
    //! \brief Function to call the doneInserting entry method of an array
    //!   fctproxy (broadcast)
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the doneInserting member function of a
    //!   chare array fctproxy and thus equivalent to
    //!   fctproxy.doneInserting(...).
    template< class Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::bcast >::value, int >::type = 0 >
    void doneDistFCTInserting( Args&&... args ) {
      fctproxy.doneInserting( std::forward<Args>(args)... );
    }

    // Calls to proxy, specific to a particular discretization

    //////  proxy.setup(...)
    //! Function to call the setup entry method of an array proxy (broadcast)
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the setup member function of a chare array
    //!   proxy and thus equivalent to proxy.setup(...), using the last argument
    //!   as default.
    template< typename... Args >
    void setup( Args&&... args ) {
      boost::apply_visitor( call_setup<Args...>( std::forward<Args>(args)... ),
                            proxy );
    }

    //////  proxy.dt(...)
    //! function to call the dt entry method of an array proxy (broadcast)
    //! \param[in] args arguments to member function (entry method) to be called
    //! \details this function calls the dt member function of a chare array
    //!   proxy and thus equivalent to proxy.dt(...), specifying a
    //!   non-default last argument.
    template< class Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::bcast >::value, int >::type = 0 >
    void dt( Args&&... args ) {
      boost::apply_visitor( call_dt<Args...>( std::forward<Args>(args)... ),
                            proxy );
    }
    //////  proxy[x].dt(...)
    //! Function to call the dt entry method of an element proxy (p2p)
    //! \param[in] x Chare array element index
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the dt member function of a chare array
    //!   element proxy and thus equivalent to proxy[x].dt(...), specifying a
    //!   non-default last argument.
    template< typename Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::elem >::value, int >::type = 0 >
    void dt( const CkArrayIndex1D& x, Args&&... args ) {
      auto e = tk::element< ProxyElem >( proxy, x );
      boost::apply_visitor( call_dt<Args...>( std::forward<Args>(args)... ),
                            e );
    }

    //////  proxy.eval(...)
    //! function to call the eval entry method of an array proxy (broadcast)
    //! \param[in] args arguments to member function (entry method) to be called
    //! \details this function calls the eval member function of a chare array
    //!   proxy and thus equivalent to proxy.eval(...), specifying a
    //!   non-default last argument.
    template< class Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::bcast >::value, int >::type = 0 >
    void eval( Args&&... args ) {
      boost::apply_visitor( call_eval<Args...>( std::forward<Args>(args)... ),
                            proxy );
    }
    //////  proxy[x].eval(...)
    //! Function to call the eval entry method of an element proxy (p2p)
    //! \param[in] x Chare array element index
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the eval member function of a chare array
    //!   element proxy and thus equivalent to proxy[x].eval(...), specifying a
    //!   non-default last argument.
    template< typename Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::elem >::value, int >::type = 0 >
    void eval( const CkArrayIndex1D& x, Args&&... args ) {
      auto e = tk::element< ProxyElem >( proxy, x );
      boost::apply_visitor( call_eval<Args...>( std::forward<Args>(args)... ),
                            e );
    }

    //////  proxy[x].insert(...)
    //! Function to call the insert entry method of an element proxy (p2p)
    //! \param[in] x Chare array element index
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the insert member function of a chare array
    //!   element proxy and thus equivalent to proxy[x].insert(...), using the
    //!   last argument as default.
    template< typename... Args >
    void insert( const CkArrayIndex1D& x, Args&&... args ) {
      auto e = tk::element< ProxyElem >( proxy, x );
      boost::apply_visitor( call_insert<Args...>( std::forward<Args>(args)... ),
                            e );
    }

    //////  proxy.doneInserting(...)
    //! \brief Function to call the doneInserting entry method of an array proxy
    //!   (broadcast)
    //! \param[in] args Arguments to member function (entry method) to be called
    //! \details This function calls the doneInserting member function of a
    //!   chare array proxy and thus equivalent to proxy.doneInserting(...).
    template< class Op, typename... Args, typename std::enable_if<
      std::is_same< Op, tag::bcast >::value, int >::type = 0 >
    void doneInserting( Args&&... args ) {
      boost::apply_visitor(
        call_doneInserting<Args...>( std::forward<Args>(args)... ), proxy );
    }

  private:
   //! Functor to call the chare entry method 'setup'
   //! \details This class is intended to be used in conjunction with variant
   //!   and boost::visitor. The template argument types are the types of the
   //!   arguments to entry method to be invoked behind the variant holding a
   //!   Charm++ proxy.
   //! \see The base class Call for the definition of operator().
   template< typename... As >
    struct call_setup : Call< call_setup<As...>, As... > {
      using Base = Call< call_setup<As...>, As... >;
      using Base::Base; // inherit base constructors
      //! Invoke the entry method
      //! \param[in,out] p Proxy behind which the entry method is called
      //! \param[in] args Function arguments passed to entry method
      //! \details P is the proxy type, Args are the types of the arguments of
      //!   the entry method to be called.
      template< typename P, typename... Args >
      static void invoke( P& p, Args&&... args ) {
        p.setup( std::forward<Args>(args)... );
      }
    };

   //! Functor to call the chare entry method 'insert'
   //! \details This class is intended to be used in conjunction with variant
   //!   and boost::visitor. The template argument types are the types of the
   //!   arguments to entry method to be invoked behind the variant holding a
   //!   Charm++ proxy.
   //! \see The base class Call for the definition of operator().
   template< typename... As >
    struct call_insert : Call< call_insert<As...>, As... > {
      using Base = Call< call_insert<As...>, As... >;
      using Base::Base; // inherit base constructors
      //! Invoke the entry method
      //! \param[in,out] p Proxy behind which the entry method is called
      //! \param[in] args Function arguments passed to entry method
      //! \details P is the proxy type, Args are the types of the arguments of
      //!   the entry method to be called.
      template< typename P, typename... Args >
      static void invoke( P& p, Args&&... args ) {
        p.insert( std::forward<Args>(args)... );
      }
    };

   //! Functor to call the chare entry method 'doneInserting'
   //! \details This class is intended to be used in conjunction with variant
   //!   and boost::visitor. The template argument types are the types of the
   //!   arguments to entry method to be invoked behind the variant holding a
   //!   Charm++ proxy.
   //! \see The base class Call for the definition of operator().
   template< typename... As >
    struct call_doneInserting : Call< call_doneInserting<As...>, As... > {
      using Base = Call< call_doneInserting<As...>, As... >;
      using Base::Base; // inherit base constructors
      //! Invoke the entry method
      //! \param[in,out] p Proxy behind which the entry method is called
      //! \param[in] args Function arguments passed to entry method
      //! \details P is the proxy type, Args are the types of the arguments of
      //!   the entry method to be called.
      template< typename P, typename... Args >
      static void invoke( P& p, Args&&... args ) {
        p.doneInserting( std::forward<Args>(args)... );
      }
    };

   //! Functor to call the chare entry method 'dt'
   //! \details This class is intended to be used in conjunction with variant
   //!   and boost::visitor. The template argument types are the types of the
   //!   arguments to entry method to be invoked behind the variant holding a
   //!   Charm++ proxy.
   //! \see The base class Call for the definition of operator().
   template< typename... As >
    struct call_dt : Call< call_dt<As...>, As... > {
      using Base = Call< call_dt<As...>, As... >;
      using Base::Base; // inherit base constructors
      //! Invoke the entry method
      //! \param[in,out] p Proxy behind which the entry method is called
      //! \param[in] args Function arguments passed to entry method
      //! \details P is the proxy type, Args are the types of the arguments of
      //!   the entry method to be called.
      template< typename P, typename... Args >
      static void invoke( P& p, Args&&... args ) {
        p.dt( std::forward<Args>(args)... );
      }
    };

   //! Functor to call the chare entry method 'eval'
   //! \details This class is intended to be used in conjunction with variant
   //!   and boost::visitor. The template argument types are the types of the
   //!   arguments to entry method to be invoked behind the variant holding a
   //!   Charm++ proxy.
   //! \see The base class Call for the definition of operator().
   template< typename... As >
    struct call_eval : Call< call_eval<As...>, As... > {
      using Base = Call< call_eval<As...>, As... >;
      using Base::Base; // inherit base constructors
      //! Invoke the entry method
      //! \param[in,out] p Proxy behind which the entry method is called
      //! \param[in] args Function arguments passed to entry method
      //! \details P is the proxy type, Args are the types of the arguments of
      //!   the entry method to be called.
      template< typename P, typename... Args >
      static void invoke( P& p, Args&&... args ) {
        p.eval( std::forward<Args>(args)... );
      }
    };

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      SchemeBase::pup( p );
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] s Scheme object reference
    friend void operator|( PUP::er& p, Scheme& s ) { s.pup(p); }
    //@}
};

} // inciter::

#endif // Scheme_h
