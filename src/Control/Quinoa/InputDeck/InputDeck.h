//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Thu 07 Aug 2014 03:28:34 PM MDT
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
                      tag::incpar,     incpars,
                      tag::component,  ncomps< int >,
                      tag::interval,   intervals,
                      tag::cmd,        CmdLine,
                      tag::param,      parameters,
                      tag::stat,       std::vector< Product > > {

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
      set< tag::incpar, tag::npar >( 1 );
      set< tag::incpar, tag::nstep >( std::numeric_limits< uint64_t >::max() );
      set< tag::incpar, tag::term >( 1.0 );
      set< tag::incpar, tag::dt >( 0.5 );
      // Default intervals
      set< tag::interval, tag::tty >( 1 );
      set< tag::interval, tag::dump >( 1 );
      set< tag::interval, tag::plot >( 1 );
      set< tag::interval, tag::pdf >( 1 );
      set< tag::interval, tag::glob >( 1 );
      // Default beta mass model parameters
      set< tag::param, tag::beta, tag::atwood >( 0.5 );
      // Default Dirichlet mix model parameters
      set< tag::param, tag::dirichlet, tag::depvar >( 'x' );
      set< tag::param, tag::dirichlet, tag::b >( std::vector< tk::real >() );
      set< tag::param, tag::dirichlet, tag::S >( std::vector< tk::real >() );
      set< tag::param, tag::dirichlet, tag::kappa >( std::vector< tk::real >() );
      // Default generalized Dirichlet mix model parameters
      set< tag::param, tag::gendir, tag::depvar >( 'x' );
      set< tag::param, tag::gendir, tag::b >( std::vector< tk::real >() );
      set< tag::param, tag::gendir, tag::S >( std::vector< tk::real >() );
      set< tag::param, tag::gendir, tag::kappa >( std::vector< tk::real >() );
      set< tag::param, tag::gendir, tag::c >( std::vector< tk::real >() );
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

  private:
    // capture map for pegtl::capture, the functions below make this class look
    // like a std::map in pegtl::capture's view, see pegtl/rules_action.hh
    pegtl::capture_map m_captured;

  public:
    //! type const_iterator for pegtl::capture
    using const_iterator = pegtl::capture_map::const_iterator;

    //! operator[] for pegtl::capture
    pegtl::capture_map::mapped_type& operator[]( unsigned key )
    { return m_captured[ key ]; }

    //! find() for pegtl::capture
    const_iterator find( const unsigned& key ) const
    { return m_captured.find( key ); }

    //! end() for pegtl::capture
    const_iterator end() const { return m_captured.end(); }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      tk::Control< tag::title,      std::string,
                   tag::selected,   selects,
                   tag::incpar,     incpars,
                   tag::component,  ncomps< int >,
                   tag::interval,   intervals,
                   tag::cmd,        CmdLine,
                   tag::param,      parameters,
                   tag::stat,       std::vector< Product > >::pup(p);
      // p | m_captured;
    }
    friend void operator|( PUP::er& p, InputDeck& c ) { c.pup(p); }
};

} // ctr::
} // quinoa::

#endif // QuinoaInputDeck_h
