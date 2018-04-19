// *****************************************************************************
/*!
  \file      src/PDE/PDEStack.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Stack of differential equations
  \details   This file declares class PDEStack, which implements various
    functionality related to registering and instantiating partial differential
    equation types. Registration and instantiation use a partial differential
    equation factory, which is a std::map (an associative container),
    associating unique partial differential equation keys to their constructor
    calls. For more details, see the in-code documentation of the constructor.
*/
// *****************************************************************************
#ifndef PDEStack_h
#define PDEStack_h

#include <map>
#include <set>
#include <string>
#include <vector>
#include <functional>

#include <boost/mpl/at.hpp>
#include <boost/mpl/aux_/adl_barrier.hpp>

#include "Tags.h"
#include "Keywords.h"
#include "Exception.h"
#include "Factory.h"
#include "CGPDE.h"
#include "DGPDE.h"
#include "Inciter/Options/PDE.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

using ncomp_t = kw::ncomp::info::expect::type;

//! \brief Factory for PDEs using continuous Galerkin discretization storing
//!   keys associated to their constructors
using CGFactory =
  std::map< ctr::PDEKey, std::function< CGPDE(const ncomp_t&) > >;

//! \brief Factory for PDEs using discontinuous Galerkin discretization storing
//!   keys associated to their constructors
using DGFactory =
  std::map< ctr::PDEKey, std::function< DGPDE(const ncomp_t&) > >;

//! \brief Partial differential equations stack
class PDEStack {

  public:
    //! Constructor: register partial differential equations into factory
    explicit PDEStack();

    //! Instantiate selected PDEs using continuous Galerkin discretization
    std::vector< CGPDE > selectedCG() const;

    //! Instantiate selected PDEs using discontinuous Galerkin discretization
    std::vector< DGPDE > selectedDG() const;

    //! Constant accessor to CGPDE factory
    //! \return Constant reference to the CGPDE factory
    const CGFactory& factory() const { return m_cgfactory; }

    //! Return info on selected partial differential equations
    std::vector< std::vector< std::pair< std::string, std::string > > > info()
    const;

    //! \brief Return number of unique equation types registered
    //! \return The number of unique equation types registered in the factory
    std::size_t ntypes() const { return m_eqTypes.size(); }

  private:
    //! \brief Function object for registering a partial differential equation
    //!   into the partial differential equation factory
    //! \details This functor is repeatedly called by MPL's cartesian_product,
    //!   sweeping all combinations of the partial differential equation
    //!   policies. The purpose of template template is to simplify client code
    //!   as that will not have to specify the template arguments of the
    //!   template argument (the policies of Eq), since we can figure it out
    //!   here. See also http://stackoverflow.com/a/214900. The template
    //!   argument Eq specifies a "child" class that is used polymorphically
    //!   with a "base" class modeling a concept defined in the base. The base
    //!   is given by the template argument PDE. The template argument Factory
    //!   specifies which factory to store the registered and configured child
    //1   PDE.
    template< template< class, class > class Eq, class Factory, class PDE >
    struct registerPDE {
      //! Need to store the reference to factory we are registering into
      Factory& factory;
      //! Need to store which differential equation we are registering
      const ctr::PDEType type;
      //! Constructor, also count number of unique equation types registered
      //! \param[in,out] f Factory into which to register PDE
      //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
      //! \param[in] eqTypes Equation type counters
      explicit registerPDE( Factory& f,
                            ctr::PDEType t,
                            std::set< ctr::PDEType >& eqTypes ) :
        factory( f ), type( t ) { eqTypes.insert( t ); }
      //! \brief Function call operator called by mpl::cartesian_product for
      //!   each unique sequence of policy combinations
      template< typename U > void operator()( U ) {
        namespace mpl = boost::mpl;
        // Get problem policy: 1st type of mpl::vector U
        using Physics = typename mpl::at< U, mpl::int_<0> >::type;
        // Get problem policy: 2nd type of mpl::vector U
        using Problem = typename mpl::at< U, mpl::int_<1> >::type;
        // Build differential equation key
        ctr::PDEKey key{ type, Physics::type(), Problem::type() };
        // Register equation (with policies given by mpl::vector U) into factory
        tk::recordModelLate< PDE, Eq< Physics, Problem > >
                           ( factory, key, static_cast<ncomp_t>(0) );
      }
    };

    //! Wrapper of registerPDE specialized for registering CG PDEs
    //! \details The sole reason for this functor is to simplify client-code
    //!   calling registerPDE specialized to CG PDEs
    template< template< class, class > class Eq >
    struct registerCG : registerPDE< Eq, CGFactory, CGPDE > {
      //! Delegate constructor to base and specialize to CG
      //! \param[in] stack Host class to access factory and counters
      //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
      explicit registerCG( PDEStack* const stack, ctr::PDEType t ) :
        registerPDE< Eq, CGFactory, CGPDE >
                   ( stack->m_cgfactory, t, stack->m_eqTypes ) {}
    };

    //! Wrapper of registerPDE specialized for registering DG PDEs
    //! \details The sole reason for this functor is to simplify client-code
    //!   calling registerPDE specialized to DG PDEs
    template< template< class, class > class Eq >
    struct registerDG : registerPDE< Eq, DGFactory, DGPDE > {
      //! Delegate constructor to base and specialize to CG
      //! \param[in] stack Host class to access factory and counters
      //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
      explicit registerDG( PDEStack* const stack, ctr::PDEType t ) :
        registerPDE< Eq, DGFactory, DGPDE >
                   ( stack->m_dgfactory, t, stack->m_eqTypes ) {}
    };

    //! \brief Instantiate a partial differential equation
    //! \details The template argument, EqTag, is used to find the given
    //!   partial differential equation in the input deck, the hierarchical data
    //!   filled during control file parsing, containing user input. The
    //!   template argument Factory specifies which factory we search in. The
    //!   template argument PDE specifies the "base" PDE type that the
    //!   instantiated "child" PDE class object is used polymorphically with.
    //! \param[in] f Factory in which to search PDE in
    //! \param[in] eq The unique partial differential equation key whose object
    //!   to instantiate.
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   partial differential equations by type.
    template< class EqTag, class Factory, class PDE >
    PDE createPDE( const Factory& f,
                   ctr::PDEType eq,
                   std::map< ctr::PDEType, ncomp_t >& cnt ) const
    {
      auto c = ++cnt[ eq ];   // count eqs
      --c;                    // used to index vectors starting with 0
      Assert( c < (g_inputdeck.get< tag::component, EqTag >().size()),
              "The number of scalar components is unspecified for the PDE to "
              "be instantiated. This is most likely a grammar error in the "
              "parser. The parser should not allow the user to select a PDE "
              "without configuring the number of scalar components the "
              "equation consists of. See inciter::deck::check_eq." );
      if ( g_inputdeck.get< tag::component, EqTag >()[c] ) {
        // re-create key and search for it
        ctr::PDEKey key{ eq,
          g_inputdeck.get< tag::param, EqTag, tag::physics >()[c],
          g_inputdeck.get< tag::param, EqTag, tag::problem >()[c] };
        const auto it = f.find( key );
        Assert( it != end( f ),
                "Can't find PDE with key('" +
                  ctr::PDE().name( key.get< tag::pde >() ) + "', '" +
                  ctr::Physics().name( key.get< tag::physics >() ) + "', '" +
                  ctr::Problem().name( key.get< tag::problem >() ) +
                  + "') in factory" );
        return it->second( c );    // instantiate and return PDE object
      } else Throw ( "Can't create PDE with zero components" );
    }

    //! Wrapper of createPDE specialized for registering CG PDEs
    //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   partial differential equations by type.
    //! \details The sole reason for this function is to simplify client-code
    //!   calling createPDE specialized to CG PDEs
    template< class EqTag >
    CGPDE createCG( ctr::PDEType t, std::map< ctr::PDEType, ncomp_t >& cnt )
    const {
      return createPDE< EqTag, CGFactory, CGPDE >( m_cgfactory, t, cnt );
    }

    //! Wrapper of createPDE specialized for registering DG PDEs
    //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   partial differential equations by type.
    //! \details The sole reason for this function is to simplify client-code
    //!   calling createPDE specialized to DG PDEs
    template< class EqTag >
    DGPDE createDG( ctr::PDEType t, std::map< ctr::PDEType, ncomp_t >& cnt )
    const {
      return createPDE< EqTag, DGFactory, DGPDE >( m_dgfactory, t, cnt );
    }

    /** @name Configuration-querying functions for PDEs */
    //! Get information on the transport PDE
    std::vector< std::pair< std::string, std::string > >
    infoTransport( std::map< ctr::PDEType, ncomp_t >& cnt ) const;
    //! Get information on the compressible flow PDEs
    std::vector< std::pair< std::string, std::string > >
    infoCompFlow( std::map< ctr::PDEType, ncomp_t >& cnt ) const;
    ///@}

    //! \brief Convert and return values from vector as string
    //! \param[in] v Vector whose components to return as a string
    //! \return Concatenated string of values read from a vector
    template< typename V >
    std::string parameters( const V& v) const {
      std::stringstream s;
      s << "{ ";
      for (auto p : v) s << p << ' ';
      s << "}";
      return s.str();
    }

    //! \brief Partial differential equations factory for those PDEs that use
    //!   continuous Galerkin discretization
    CGFactory m_cgfactory;
    //! \brief Partial differential equations factory for those PDEs that use
    //!   discontinuous Galerkin discretization
    DGFactory m_dgfactory;
    //! Counters for equation types
    std::set< ctr::PDEType > m_eqTypes;
};

} // inciter::

#endif // PDEStack_h
