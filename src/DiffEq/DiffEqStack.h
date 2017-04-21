// *****************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Stack of differential equations
  \details   This file declares class DiffEqStack, which implements various
    functionality related to registering and instantiating differential equation
    types. Registration and instantiation use a differential equation factory,
    which is a std::map (an associative container), associating unique
    differential equation keys to their constructor calls. For more details, see
    the in-code documentation of the constructor.
*/
// *****************************************************************************
#ifndef DiffEqStack_h
#define DiffEqStack_h

#include <map>
#include <set>
#include <string>
#include <vector>
#include <functional>
#include <iterator>
#include <ostream>
#include <utility>
#include <cstddef>

#include <boost/mpl/at.hpp>
#include <boost/mpl/aux_/adl_barrier.hpp>

#include "Tags.h"
#include "Keywords.h"
#include "Exception.h"
#include "Factory.h"
#include "DiffEq.h"
#include "SystemComponents.h"
#include "Walker/Options/DiffEq.h"
#include "Walker/InputDeck/InputDeck.h"

namespace walker {

extern ctr::InputDeck g_inputdeck;

//! \brief Differential equation factory: keys associated to their constructors
//! \author J. Bakosi
using DiffEqFactory =
  std::map< ctr::DiffEqKey,
            std::function< DiffEq(const tk::ctr::ncomp_type&) > >;

//! \brief Differential equations stack
//! \author J. Bakosi
class DiffEqStack {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! Constructor: register differential equations into factory
    explicit DiffEqStack();

    //! Instantiate selected DiffEqs
    std::vector< DiffEq > selected() const;

    //! \brief Instantiate tables from which extra statistics data to be output
    //!    sampled for all selected differential equations
    std::pair< std::vector< std::string >, std::vector< tk::Table > > tables()
    const;

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
                           ( factory, key, static_cast<ncomp_t>(0) );
      }
    };

    //! \brief Instantiate a differential equation object
    //! \see walker::DiffEq
    //! \details Since multiple systems of equations can be configured
    //!   with the same type, the map (cnt) is used to assign a counter to each
    //!   type. The template argument, EqTag, is used to find the given system
    //!   of differential equations in the input deck, the hierarchical data
    //!   filled during control file parsing, containing user input.
    //! \param[in] eq The unique differential equation system key whose object
    //!   to instantiate.
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   differential equations systems by type.
    //! \author J. Bakosi
    template< class EqTag >
    DiffEq createDiffEq( ctr::DiffEqType eq,
                         std::map< ctr::DiffEqType, ncomp_t >& cnt ) const {
      auto c = ++cnt[ eq ];   // count eqs
      --c;                    // used to index vectors starting with 0
      if ( g_inputdeck.get< tag::component, EqTag >()[c] ) {
        // create key and search for it
        ctr::DiffEqKey key{ eq,
          g_inputdeck.get< tag::param, EqTag, tag::initpolicy >()[c],
          g_inputdeck.get< tag::param, EqTag, tag::coeffpolicy >()[c] };
        const auto it = m_factory.find( key );
        Assert( it != end( m_factory ), "Can't find eq in factory" );
        // instantiate and return diff eq object
        return it->second( c );
      } else Throw ( "Can't create DiffEq with zero independent variables" );
    }

    //! \brief Instantiate tables from a differential equation
    //! \details This function is used to instantiate a vector of tk::Tables and
    //!   their associated names for a system of differential equation selected
    //!   by the user. Since multiple systems of equations can be configured
    //!   with the same type, the map (cnt) is used to assign a counter to each
    //!   type. We then find out if the coefficients policy, configured to the
    //!   equation system, uses tables. If so, we return all the tables and
    //!   their names associated to the different components in the system. The
    //!   template argument, EqTag, is used to find the given differential
    //!   equation system in the input deck, the hierarchical data filled
    //!   during control file parsing, containing user input.
    //! \param[in] eq The unique differential equation system key whose tables
    //!   object to instantiate (if the configured coefficients policy uses
    //!   tables).
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   vector of tables associated to a differential equations systems by
    //!   type.
    //! \author J. Bakosi
    template< class EqTag >
    std::pair< std::vector< std::string >, std::vector< tk::Table > >
    createTables( ctr::DiffEqType eq,
                  std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
    {
      auto c = ++cnt[ eq ];   // count eqs
      --c;                    // used to index vectors starting with 0
      std::vector< tk::Table > tab;
      std::vector< std::string > nam;
      if ( g_inputdeck.get< tag::component, EqTag >()[c] ) {
        // find out if coefficients policy uses tables and return them if so
        if (g_inputdeck.get< tag::param, EqTag, tag::coeffpolicy >()[c] ==
              ctr::CoeffPolicyType::HYDROTIMESCALE_HOMOGENEOUS_DECAY)
        {
          const auto& hts = g_inputdeck.get< tag::param,
                                             tag::mixmassfracbeta,
                                             tag::hydrotimescales >().at(c);
          ctr::HydroTimeScales ot;
          for (auto t : hts) tab.push_back( ot.table(t) );
          for (auto t : hts) nam.push_back( ot.name(t) );

          const auto& hp = g_inputdeck.get< tag::param,
                                            tag::mixmassfracbeta,
                                            tag::hydroproductions >().at(c);
          ctr::HydroProductions op;
          for (auto t : hp) tab.push_back( op.table(t) );
          for (auto t : hp) nam.push_back( op.name(t) );

        }
      } else Throw ( "DiffEq with zero independent variables" );
      return { nam, tab };
    }

    /** @name Configuration-querying functions for SDEs */
    //! Get information on the Dirichlet SDE
    std::vector< std::pair< std::string, std::string > >
    infoDirichlet( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on Lochner's generalized Dirichlet SDE
    std::vector< std::pair< std::string, std::string > >
    infoGenDir( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on Wright-Fisher SDE
    std::vector< std::pair< std::string, std::string > >
    infoWrightFisher( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on Ornstein_Uhlenbeck SDE
    std::vector< std::pair< std::string, std::string > >
    infoOU( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on diagonal Ornstein_Uhlenbeck SDE
    std::vector< std::pair< std::string, std::string > >
    infoDiagOU( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on beta SDE
    std::vector< std::pair< std::string, std::string > >
    infoBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on number-fraction beta SDE
    std::vector< std::pair< std::string, std::string > >
    infoNumberFractionBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on mass-fraction beta SDE
    std::vector< std::pair< std::string, std::string > >
    infoMassFractionBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on mix number-fraction beta SDE
    std::vector< std::pair< std::string, std::string > >
    infoMixNumFracBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on mix mass-fraction beta SDE
    std::vector< std::pair< std::string, std::string > >
    infoMixMassFracBeta( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on skew-normal SDE
    std::vector< std::pair< std::string, std::string > >
    infoSkewNormal( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
    //! Get information on Gamma SDE
    std::vector< std::pair< std::string, std::string > >
    infoGamma( std::map< ctr::DiffEqType, ncomp_t >& cnt ) const;
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

    //! \brief Return names of options (tk::Toggle) from vector as a string
    //! \param[in] v Option vector whose names of components to return
    //! \return Concatenated string of option names read from option vector
    //! \author J. Bakosi
    template< class Option, class OptTypeVec >
    std::string options( const Option& opt, const OptTypeVec& v ) const {
      std::stringstream s;
      s << "{ ";
      for (auto o : v) s << opt.name(o) << ' ';
      s << "}";
      return s.str();
    }

    //! \brief Insert spike information (used to specify delta PDFs) into info
    //!   vector
    //! \param[in,out] nfo Info vector of string-pairs to insert to
    //! \param[in] spike Vector of vectors specifying spike info
    //! \author J. Bakosi
    template< typename Info, typename VV >
    void spikes( Info& nfo, const VV& spike ) const {
      std::size_t i = 0;
      for (const auto& s : spike)
        nfo.emplace_back( "delta spikes [comp" + std::to_string(++i) + ":" +
                            std::to_string( s.size()/2 ) + "]",
                          parameters( s ) );
    }

    //! \brief Insert betapdf information (used to specify beta PDFs) into info
    //!   vector
    //! \param[in,out] nfo Info vector of string-pairs to insert to
    //! \param[in] betapdf Vector of vectors specifying betapdf info
    //! \author J. Bakosi
    template< typename Info, typename VV >
    void betapdfs( Info& nfo, const VV& betapdf ) const {
      std::size_t i = 0;
      for (const auto& s : betapdf)
        nfo.emplace_back( "beta pds [comp" + std::to_string(++i) + "]",
                          parameters( s ) );
    }

    DiffEqFactory m_factory;                 //!< Differential equations factory
    std::set< ctr::DiffEqType > m_eqTypes;   //!< Count number of equation types
};

} // walker::

#endif // DiffEqStack_h
