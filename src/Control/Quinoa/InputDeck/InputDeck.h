//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Thu 04 Sep 2014 04:57:25 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's input deck
  \details   Quinoa's input deck
*/
//******************************************************************************
#ifndef QuinoaInputDeck_h
#define QuinoaInputDeck_h

#include <limits>

#include <Control.h>
#include <Quinoa/CmdLine/CmdLine.h>
#include <Quinoa/Components.h>

namespace quinoa {
namespace ctr {

//! InputDeck : Control< specialized to Quinoa >, see Types.h,
class InputDeck :
  public tk::Control< // tag           type
                      tag::title,      std::string,
                      tag::selected,   selects,
                      tag::discr,      discretization,
                      tag::component,  ncomps< unsigned int >,
                      tag::interval,   intervals,
                      tag::cmd,        CmdLine,
                      tag::param,      parameters,
                      tag::stat,       std::vector< Product >,
                      tk::tag::error,  std::vector< std::string > > {

  public:
    //! Constructor: set all defaults
    InputDeck() {
      // Default title
      set< tag::title >( "" );
      // Default selections
      set< tag::selected, tag::montecarlo >( MonteCarloType::NO_MONTECARLO );
      set< tag::selected, tag::position >( PositionType::NO_POSITION );
      set< tag::selected, tag::mass >( MassType::NO_MASS );
      set< tag::selected, tag::hydro >( HydroType::NO_HYDRO );
      set< tag::selected, tag::energy >( EnergyType::NO_ENERGY );
      set< tag::selected, tag::mix >( MixType::NO_MIX );
      set< tag::selected, tag::frequency >( FrequencyType::NO_FREQUENCY );
      set< tag::selected, tag::mixrate >( MixRateType::NO_MIXRATE );
      // Default time incrementation parameters
      set< tag::discr, tag::npar >( 1 );
      set< tag::discr, tag::nstep >( std::numeric_limits< uint64_t >::max() );
      set< tag::discr, tag::term >( 1.0 );
      set< tag::discr, tag::dt >( 0.5 );
      set< tag::discr, tag::dt >( 0.5 );
      // Default intervals
      set< tag::interval, tag::tty >( 1 );
      set< tag::interval, tag::dump >( 1 );
      set< tag::interval, tag::plot >( 1 );
      set< tag::interval, tag::pdf >( 1 );
      set< tag::interval, tag::glob >( 1 );
      // Default beta mass model parameters
      set< tag::param, tag::beta, tag::atwood >( 0.5 );
      // Default gamma mix model parameters
      set< tag::param, tag::gamma, tag::c1 >( 0.5 );
      set< tag::param, tag::gamma, tag::c2 >( 0.73 );
      set< tag::param, tag::gamma, tag::c3 >( 5.0 );
      set< tag::param, tag::gamma, tag::c4 >( 0.25 );
      // Default simplified Langevin hydro model parameters
      set< tag::param, tag::slm, tag::c0 >( 2.1 );
      // Default generalized Langevin hydro model parameters
      set< tag::param, tag::slm, tag::c0 >( 2.1 );
      // Default requested statistics
      set< tag::stat >( std::vector< Product >() );
    }

    //! Extract data on whether to plot ordinary moments of requested statistics
    std::vector< bool > plotOrdinary() const {
      std::vector< bool > plot;
      for (const auto& product : get< tag::stat >()) {
        if (ordinary( product )) {
          plot.emplace_back( false );
          for (const auto& term : product)
            if (term.plot) plot.back() = true;
        }
      }
      return plot;
    }

    //! Extract moment names of requested statistics
    std::vector< std::string > momentNames( std::function<
      bool ( const std::vector< ctr::Term >& ) > momentType ) const
    {
      std::vector< std::string > names;
      for (const auto& product : get< tag::stat >()) {
        if (momentType( product )) {
          names.emplace_back( std::string() );
          for (const auto& term : product)
            names.back() += ctr::FieldVar( term.var, term.field );
        }
      }
      return names;
    }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      tk::Control< tag::title,      std::string,
                   tag::selected,   selects,
                   tag::discr,      discretization,
                   tag::component,  ncomps< unsigned int >,
                   tag::interval,   intervals,
                   tag::cmd,        CmdLine,
                   tag::param,      parameters,
                   tag::stat,       std::vector< Product >,
                   tk::tag::error,  std::vector< std::string > >::pup(p);
    }
    friend void operator|( PUP::er& p, InputDeck& c ) { c.pup(p); }
};

} // ctr::
} // quinoa::

#endif // QuinoaInputDeck_h
