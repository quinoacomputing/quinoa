// *****************************************************************************
/*!
  \file      src/DiffEq/DiffEqStack.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "NoWarning/back.hpp"
#include "NoWarning/front.hpp"

#include "Tags.hpp"
#include "Exception.hpp"
#include "DiffEq.hpp"
#include "DiffEqFactory.hpp"
#include "SystemComponents.hpp"
#include "Walker/InputDeck/InputDeck.hpp"

namespace walker {

extern ctr::InputDeck g_inputdeck;

//! \brief Differential equations stack
class DiffEqStack {

  private:
    using ncomp_t = tk::ctr::ncomp_t;

  public:
    //! Constructor: register differential equations into factory
    explicit DiffEqStack();

    //! Instantiate selected DiffEqs
    std::vector< DiffEq > selected() const;

    //! \brief Instantiate tables from which extra statistics data to be output
    //!    sampled for all selected differential equations
    std::pair< std::vector< std::string >, std::vector< tk::Table<1> > >
    tables() const;

    //! \brief Constant accessor to differential equation factory
    //! \return Constant reference to the internal differential equation factory
    const DiffEqFactory& factory() const { return m_factory; }

    //! Return info on selected differential equations
    std::vector< std::vector< std::pair< std::string, std::string > > > info()
    const;

    //! \brief Return number of unique equation types registered
    //! \return The number of unique equation types registered in the factory
    std::size_t ntypes() const { return m_eqTypes.size(); }

  private:
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
    template< class EqTag >
    DiffEq createDiffEq( ctr::DiffEqType eq,
                         std::map< ctr::DiffEqType, ncomp_t >& cnt ) const {
      auto c = ++cnt[ eq ];   // count eqs
      --c;                    // used to index vectors starting with 0
      if ( g_inputdeck.get< tag::component, EqTag >()[c] ) {
        // create key and search for it
        ctr::DiffEqKey key{{ eq,
          g_inputdeck.get< tag::param, EqTag, tag::initpolicy >()[c],
          g_inputdeck.get< tag::param, EqTag, tag::coeffpolicy >()[c] }};
        const auto it = m_factory.find( key );
        Assert( it != end( m_factory ),
                "Can't find eq '" + ctr::DiffEq().name( eq ) +
                "' in DiffEq factory with initialization policy '" +
                ctr::InitPolicy().name(
                  g_inputdeck.get< tag::param, EqTag, tag::initpolicy >()[c] ) +
                "' and coefficient policy '" +
                ctr::CoeffPolicy().name(
                  g_inputdeck.get< tag::param, EqTag, tag::coeffpolicy >()[c] )
                + "'" );
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
    template< class EqTag >
    std::pair< std::vector< std::string >, std::vector< tk::Table<1> > >
    createTables( ctr::DiffEqType eq,
                  std::map< ctr::DiffEqType, ncomp_t >& cnt ) const
    {
      auto c = ++cnt[ eq ];   // count eqs
      --c;                    // used to index vectors starting with 0
      std::vector< tk::Table<1> > tab;
      std::vector< std::string > nam;
      const auto& ncompeq = g_inputdeck.get< tag::component, EqTag >();
      if (!ncompeq.empty()) {
        if ( g_inputdeck.get< tag::component, EqTag >()[c] ) {
          // find out if coefficients policy uses tables and return them if so
          if (g_inputdeck.get< tag::param, EqTag, tag::coeffpolicy >()[c] ==
                ctr::CoeffPolicyType::HYDROTIMESCALE)
          {
            const auto& hts = g_inputdeck.get< tag::param,
                                               EqTag,
                                               tag::hydrotimescales >().at(c);
            ctr::HydroTimeScales ot;
            // cppcheck-suppress useStlAlgorithm
            for (auto t : hts) tab.push_back( ot.table(t) );
            // cppcheck-suppress useStlAlgorithm
            for (auto t : hts) nam.push_back( ot.name(t) );
  
            const auto& hp = g_inputdeck.get< tag::param,
                                              EqTag,
                                              tag::hydroproductions >().at(c);
            ctr::HydroProductions op;
            // cppcheck-suppress useStlAlgorithm
            for (auto t : hp) tab.push_back( op.table(t) );
            // cppcheck-suppress useStlAlgorithm
            for (auto t : hp) nam.push_back( op.name(t) );
  
          }
        } else Throw ( "DiffEq with zero independent variables" );
      }
      return { nam, tab };
    }

    DiffEqFactory m_factory;                 //!< Differential equations factory
    std::set< ctr::DiffEqType > m_eqTypes;   //!< Count number of equation types
};

} // walker::

#endif // DiffEqStack_h
