// *****************************************************************************
/*!
  \file      src/PDE/PDEStack.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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
#include "PDE.h"
#include "Inciter/Options/PDE.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

using ncomp_t = kw::ncomp::info::expect::type;

//! \brief Partial differential equation factory: keys associated to their
//!   constructors
//! \author J. Bakosi
using PDEFactory =
  std::map< ctr::PDEKey, std::function< PDE(const ncomp_t&) > >;

//! \brief Partial differential equations stack
//! \author J. Bakosi
class PDEStack {

  public:
    //! Constructor: register partial differential equations into factory
    explicit PDEStack();

    //! Instantiate selected partial differential equations
    std::vector< PDE > selected() const;

    //! \brief Constant accessor to partial differential equation factory
    //! \return Constant reference to the internal partial differential equation
    //!   factory
    //! \author J. Bakosi
    const PDEFactory& factory() const { return m_factory; }

    //! Return info on selected partial differential equations
    std::vector< std::vector< std::pair< std::string, std::string > > > info()
    const;

    //! \brief Return number of unique equation types registered
    //! \return The number of unique equation types registered in the factory
    //! \author J. Bakosi
    std::size_t ntypes() const { return m_eqTypes.size(); }

  private:
    //! \brief Function object for registering a partial differential equation
    //!   into the partial differential equation factory
    //! \details This functor is repeatedly called by MPL's cartesian_product,
    //!   sweeping all combinations of the partial differential equation
    //!   policies. The purpose of template template is to simplify client code
    //!   as that will not have to specify the template arguments of the
    //!   template argument (the policies of Eq), since we can figure it out
    //!   here. See also http://stackoverflow.com/a/214900
    //! \author J. Bakosi
    template< template< class, class > class Eq >
    struct registerPDE {
      //! Need to store the reference to factory we are registering into
      PDEFactory& factory;
      //! Need to store which differential equation we are registering
      const ctr::PDEType type;
      //! Constructor, also count number of unique equation types registered
      explicit registerPDE( PDEFactory& f,
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

    //! \brief Instantiate a partial differential equation
    //! \details The template argument, EqTag, is used to find the given
    //!   partial differential equation in the input deck, the hierarchical data
    //!   filled during control file parsing, containing user input.
    //! \param[in] eq The unique partial differential equation key whose object
    //!   to instantiate.
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   partial differential equations by type.
    //! \author J. Bakosi
    template< class EqTag >
    PDE createPDE( ctr::PDEType eq, std::map< ctr::PDEType, ncomp_t >& cnt )
    const {
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
        const auto it = m_factory.find( key );
        Assert( it != end( m_factory ), "Can't find PDE in factory" );
        return it->second( c );    // instantiate and return PDE object
      } else Throw ( "Can't create PDE with zero components" );
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
    //! \author J. Bakosi
    template< typename V >
    std::string parameters( const V& v) const {
      std::stringstream s;
      s << "{ ";
      for (auto p : v) s << p << ' ';
      s << "}";
      return s.str();
    }

    PDEFactory m_factory;            //!< Partial differential equations factory
    std::set< ctr::PDEType > m_eqTypes;   //!< Count number of equation types
};

} // inciter::

#endif // PDEStack_h
