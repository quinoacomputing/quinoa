//******************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.h
  \author    J. Bakosi
  \date      Wed 17 Dec 2014 04:41:14 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Stack of differential equations
  \details   Stack of differential equations
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

//! Differential equation factory: keys associated to their constructors
using DiffEqFactory = std::map< ctr::DiffEqKey, std::function< DiffEq(int) > >;

//! DiffEqStack
class DiffEqStack {

  public:
    //! Constructor: register differential equations into factory
    explicit DiffEqStack();

    //! Instantiate selected DiffEqs
    std::vector< DiffEq > selected() const;

    //! Accessor to factory
    const DiffEqFactory& factory() const { return m_factory; }

    //! Return info on selected differential equations
    std::vector< std::vector< std::pair< std::string, std::string > > > info()
    const;

    //! Return number of unique equation types registered
    std::size_t ntypes() const { return m_eqTypes.size(); }

  private:
    //! Function object for registering a differential equation into the
    //! differential equation factory - repeatedly called by mpl's
    //! cartesian_product sweeping all combinations of the differential equation
    //! policies. The purpose of template template is to simplify client code as
    //! that will not have to specify the template arguments of the template
    //! argument (the policies of Eq), since we can figure it out here. See also
    //! http://stackoverflow.com/a/214900
    template< template< class, class > class Eq >
    struct registerDiffEq {

      DiffEqFactory& factory;
      const ctr::DiffEqType type;

      // Constructor, also count number of unique equation types registered
      explicit registerDiffEq( DiffEqFactory& f, ctr::DiffEqType t,
                               std::set< ctr::DiffEqType >& eqTypes ) :
        factory( f ), type( t ) { eqTypes.insert( t ); }

      template< typename U > void operator()( U ) {
        namespace mpl = boost::mpl;
        // Get Initialization policy: 1st type of mpl::vector U
        using InitPolicy = typename mpl::at< U, mpl::int_<0> >::type;
        // Get coefficients policy: 2nd type of mpl::vector U
        using CoeffPolicy = typename mpl::at< U, mpl::int_<1> >::type;
        // Build differential equation key
        ctr::DiffEqKey key{ type, InitPolicy().type(), CoeffPolicy().type() };
        // Register equation (with policies given by mpl::vector U) into factory
        tk::recordModelLate< DiffEq, Eq< InitPolicy, CoeffPolicy > >
                           ( factory, key, 0 );
      }
    };

    //! Instantiate differential equation
    template< class EqTag >
    DiffEq createDiffEq( const DiffEqFactory& f,
                         ctr::DiffEqType eq,
                         std::map< ctr::DiffEqType, int >& cnt ) const {
      auto c = ++cnt[ eq ];   // count eqs
      --c;                    // used to index vectors starting with 0
      if ( g_inputdeck.get< tag::component, EqTag >()[c] ) {
        ctr::DiffEqKey key{ eq,
          g_inputdeck.get< tag::param, EqTag, tag::initpolicy >()[c],
          g_inputdeck.get< tag::param, EqTag, tag::coeffpolicy >()[c] };
        const auto it = f.find( key );
        Assert( it != end( f ), "Can't find eq in factory" );
        return it->second( c );    // instantiate and return diff eq
      } else Throw ( "Can't create DiffEq with zero independent variables" );
    }

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

    //! Return parameter values from vector as string
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
