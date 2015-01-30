//******************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 05:49:44 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Stack of differential equations
  \details   This file declares class DiffEqStack, which implements various
    functionality related to registering and instantiating differential equation
    types. Registration and instantiation use a differential equation factory,
    which is a std::map (an associative container), associating unique
    differential equation keys to their constructor calls. For more details, see
    the in-code documentation of the constructor.
*/
//******************************************************************************
#ifndef DiffEqStack_h
#define DiffEqStack_h

#include <vector>

#include <boost/mpl/at.hpp>

#include <DiffEq.h>
#include <Factory.h>
#include <Tags.h>
#include <Walker/Options/DiffEq.h>
#include <Options/InitPolicy.h>
#include <Options/CoeffPolicy.h>
#include <Walker/InputDeck/InputDeck.h>

namespace walker {

extern ctr::InputDeck g_inputdeck;

//! \brief Differential equation factory: keys associated to their constructors
//! \author J. Bakosi
using DiffEqFactory = std::map< ctr::DiffEqKey, std::function< DiffEq(int) > >;

//! \brief Differential equations stack
//! \author J. Bakosi
class DiffEqStack {

  public:
    //! Constructor: register differential equations into factory
    explicit DiffEqStack();

    //! Instantiate selected DiffEqs
    std::vector< DiffEq > selected() const;

    //! \brief Constant accessor to differential equation factory
    //! \return Constant reference to the internal differential equation factory
    //! \author J. Bakosi
    const DiffEqFactory& factory() const { return m_factory; }

    //! Return info on selected differential equations
    std::vector< std::vector< std::pair< std::string, std::string > > > info()
    const;

    //! \brief Return number of unique equation types registered
    //! \return The number of unique equation types registered in the factory
    //! \author J. Bakosi
    std::size_t ntypes() const { return m_eqTypes.size(); }

  private:
    //! \brief Function object for registering a differential equation into the
    //!   differential equation factory
    //! \details This functor is repeatedly called by MPL's cartesian_product,
    //!   sweeping all combinations of the differential equation policies. The
    //!   purpose of template template is to simplify client code as that will
    //!   not have to specify the template arguments of the template argument
    //!   (the policies of Eq), since we can figure it out here. See also
    //!   http://stackoverflow.com/a/214900
    //! \author J. Bakosi
    template< template< class, class > class Eq >
    struct registerDiffEq {
      //! Need to store the reference to factory we are registering into
      DiffEqFactory& factory;
      //! Need to store which differential equation we are registering
      const ctr::DiffEqType type;
      //! Constructor, also count number of unique equation types registered
      explicit registerDiffEq( DiffEqFactory& f,
                               ctr::DiffEqType t,
                               std::set< ctr::DiffEqType >& eqTypes ) :
        factory( f ), type( t ) { eqTypes.insert( t ); }
      //! \brief Function call operator called by mpl::cartesian_product for
      //!   each unique sequence of policy combinations
      template< typename U > void operator()( U ) {
        namespace mpl = boost::mpl;
        // Get Initialization policy: 1st type of mpl::vector U
        using InitPolicy = typename mpl::at< U, mpl::int_<0> >::type;
        // Get coefficients policy: 2nd type of mpl::vector U
        using CoeffPolicy = typename mpl::at< U, mpl::int_<1> >::type;
        // Build differential equation key
        ctr::DiffEqKey key{ type, InitPolicy::type(), CoeffPolicy::type() };
        // Register equation (with policies given by mpl::vector U) into factory
        tk::recordModelLate< DiffEq, Eq< InitPolicy, CoeffPolicy > >
                           ( factory, key, 0 );
      }
    };

    //! \brief Instantiate a differential equation
    //! \details The template argument, EqTag, is used to find the given
    //!   differential equation in the input deck, the hierarchical data filled
    //!   during control file parsing, containing user input.
    //! \param[in] eq The unique differential equation key whose object to
    //!   instantiate.
    //! \param[inout] cnt Counter, a std::map, that counts all instantiated
    //!   differential equations by type.
    //! \author J. Bakosi
    template< class EqTag >
    DiffEq createDiffEq( ctr::DiffEqType eq,
                         std::map< ctr::DiffEqType, int >& cnt ) const {
      auto c = ++cnt[ eq ];   // count eqs
      --c;                    // used to index vectors starting with 0
      if ( g_inputdeck.get< tag::component, EqTag >()[c] ) {
        // re-create key and search for it
        ctr::DiffEqKey key{ eq,
          g_inputdeck.get< tag::param, EqTag, tag::initpolicy >()[c],
          g_inputdeck.get< tag::param, EqTag, tag::coeffpolicy >()[c] };
        const auto it = m_factory.find( key );
        Assert( it != end( m_factory ), "Can't find eq in factory" );
        return it->second( c );    // instantiate and return diff eq object
      } else Throw ( "Can't create DiffEq with zero independent variables" );
    }

    /** @name Configuration-querying functions for all SDEs */
    //! Get information on the Dirichlet SDE
    std::vector< std::pair< std::string, std::string > >
    infoDirichlet( std::map< ctr::DiffEqType, int >& cnt ) const;
    //! Get information on Lochner's generalized Dirichlet SDE
    std::vector< std::pair< std::string, std::string > >
    infoGenDir( std::map< ctr::DiffEqType, int >& cnt ) const;
    //! Get information on Wright-Fisher SDE
    std::vector< std::pair< std::string, std::string > >
    infoWrightFisher( std::map< ctr::DiffEqType, int >& cnt ) const;
    //! Get information on Ornstein_Uhlenbeck SDE
    std::vector< std::pair< std::string, std::string > >
    infoOU( std::map< ctr::DiffEqType, int >& cnt ) const;
    //! Get information on diagonal Ornstein_Uhlenbeck SDE
    std::vector< std::pair< std::string, std::string > >
    infoDiagOU( std::map< ctr::DiffEqType, int >& cnt ) const;
    //! Get information on Beta SDE
    std::vector< std::pair< std::string, std::string > >
    infoBeta( std::map< ctr::DiffEqType, int >& cnt ) const;
    //! Get information on skew-normal SDE
    std::vector< std::pair< std::string, std::string > >
    infoSkewNormal( std::map< ctr::DiffEqType, int >& cnt ) const;
    //! Get information on Gamma SDE
    std::vector< std::pair< std::string, std::string > >
    infoGamma( std::map< ctr::DiffEqType, int >& cnt ) const;
    ///@}

    //! \brief Return parameter values from vector as string
    //! \details The template arguments index into the input deck, the
    //!   hierarchical data structure filled during control file parsing,
    //!   containing user input.
    //! \param[in] c Index into a vector to access a vector
    //! \return Concatenated string of values read from a vector
    //! \author J. Bakosi
    template< typename... tags > std::string parameters( int c ) const {
      std::stringstream s;
      s << "{ ";
      for (auto p : g_inputdeck.get< tags... >()[c]) s << p << ' ';
      s << "}";
      return s.str();
    }

    DiffEqFactory m_factory;                 //!< Differential equations factory
    std::set< ctr::DiffEqType > m_eqTypes;   //!< Count number of equation types
};

} // walker::

#endif // DiffEqStack_h
