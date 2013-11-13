//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Tue 12 Nov 2013 09:24:35 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck
  \details   Quinoa's input deck
*/
//******************************************************************************
#ifndef QuinoaInputDeck_h
#define QuinoaInputDeck_h

#include <limits>

#include <Control.h>
#include <Option.h>
#include <Quinoa/CmdLine/CmdLine.h>

namespace quinoa {
namespace ctr {

//! InputDeck : Control< specialized to Quinoa >, see Types.h,
class InputDeck :
  public tk::Control< // tag      type
                      title,      std::string,
                      selected,   selects,
                      incpar,     incpars,
                      component,  components,
                      interval,   intervals,
                      cmd,        CmdLine,
                      param,      parameters,
                      stat,       std::vector< Product > > {

  public:
    //! Constructor: set all defaults
    InputDeck() {
      // Default title
      set< ctr::title >( "" );
      // Default selections
      set< selected, geometry >( GeometryType::NO_GEOMETRY );
      set< selected, physics >( PhysicsType::NO_PHYSICS );
      set< selected, position >( PositionType::NO_POSITION );
      set< selected, mass >( MassType::NO_MASS );
      set< selected, hydro >( HydroType::NO_HYDRO );
      set< selected, energy >( EnergyType::NO_ENERGY );
      set< selected, mix >( MixType::NO_MIX );
      set< selected, frequency >( FrequencyType::NO_FREQUENCY );
      set< selected, mixrate >( MixRateType::NO_MIXRATE );
      set< selected, rng >( RNGType::NO_RNG );
      // Default time incrementation parameters
      set< incpar, nstep >( std::numeric_limits< uint64_t >::max() );
      set< incpar, term >( 1.0 );
      set< incpar, dt >( 0.5 );
      // Default number of components
      set< component, nposition >( 0 );
      set< component, ndensity >( 0 );
      set< component, nvelocity >( 0 );
      set< component, nscalar >( 0 );
      set< component, nfrequency >( 0 );
      set< component, npar >( 1 );
      // Default intervals
      set< interval, tty >( 1 );
      set< interval, dump >( 1 );
      set< interval, plot >( 1 );
      set< interval, pdf >( 1 );
      set< interval, glob >( 1 );
      // Default random number generator seed
      set< param, rng, seed >( 0 );
      // Default beta mass model parameters
      set< param, beta, atwood >( 0.5 );
      // Default Dirichlet mix model parameters
      set< param, dirichlet, b >( std::vector< tk::real >() );
      set< param, dirichlet, S >( std::vector< tk::real >() );
      set< param, dirichlet, kappa >( std::vector< tk::real >() );
      // Default generalized Dirichlet mix model parameters
      set< param, gendirichlet, b >( std::vector< tk::real >() );
      set< param, gendirichlet, S >( std::vector< tk::real >() );
      set< param, gendirichlet, kappa >( std::vector< tk::real >() );
      set< param, gendirichlet, c >( std::vector< tk::real >() );
      // Default gamma mix model parameters
      set< param, gamma, c1 >( 0.5 );
      set< param, gamma, c2 >( 0.73 );
      set< param, gamma, c3 >( 5.0 );
      set< param, gamma, c4 >( 0.25 );
      // Default simplified Langevin hydro model parameters
      set< param, slm, c0 >( 2.1 );
      // Default generalized Langevin hydro model parameters
      set< param, slm, c0 >( 2.1 );
      // Default requested statistics
      set<stat>( std::vector< Product >() );
    }

    //! Return total number of particle properties
    uint32_t nprop() const noexcept {
      return get< component, nposition >() +
             get< component, ndensity >() +
             get< component, nvelocity> () +
             get< component, nscalar >() +
             get< component, nfrequency >();
    }

    //! Return position offset
    int positionOffset() const noexcept {
      return 0;
    }
    //! Return density offset
    int densityOffset() const noexcept {
      return get< component, nposition >();
    }
    //! Return velocity offset
    int velocityOffset() const noexcept {
      return get< component, nposition >() +
             get< component, ndensity >();
    }
    //! Return scalar offset
    int scalarOffset() const noexcept {
      return get< component, nposition >() +
             get< component, ndensity >() +
             get< component, nvelocity >();
    }

    //! Return offset for term::quantity
    int termOffset(Quantity q) const noexcept {
      int offset = 0;
      if (q == Quantity::SCALAR)
        offset += get< component, nvelocity >();
      if (q == Quantity::VELOCITY_Z)
        offset += get< component, nvelocity >();
      if (q == Quantity::VELOCITY_Y)
        offset += get< component, nvelocity >();
      if (q == Quantity::VELOCITY_X)
        offset += get< component, ndensity >();
      if (q == Quantity::DENSITY)
        offset += 3 * get< component, nposition >();
      return offset;
    }

  private:
    //! Don't permit copy constructor
    InputDeck(const InputDeck&) = delete;
    //! Don't permit copy assigment
    InputDeck& operator=(const InputDeck&) = delete;
    //! Don't permit move constructor
    InputDeck(InputDeck&&) = delete;
    //! Don't permit move assigment
    InputDeck& operator=(InputDeck&&) = delete;
};

//! InputDeck defaults
static const InputDeck InputDeckDefaults;

} // ctr::
} // quinoa::

#endif // QuinoaInputDeck_h
