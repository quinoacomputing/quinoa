//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Thu Oct  3 17:13:06 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck
  \details   Quinoa's input deck
*/
//******************************************************************************
#ifndef QuinoaInputDeck_h
#define QuinoaInputDeck_h

#include <Control.h>
#include <Option.h>
#include <Quinoa/InputDeck/Types.h>

namespace quinoa {

//! InputDeck : Control< specialized to Quinoa >, see Types.h,
//! This is also where the command line parser stores
class InputDeck :
  public Control< // tag           type
                  ctr::title,      std::string,
                  ctr::selected,   ctr::selects,
                  ctr::incpar,     ctr::incpars,
                  ctr::component,  ctr::components,
                  ctr::interval,   ctr::intervals,
                  ctr::io,         ctr::ios,
                  ctr::param,      ctr::parameters,
                  ctr::stat,       std::vector<ctr::Product> > {

  public:
    //! Constructor: set all defaults
    InputDeck() {
      using namespace ctr;
      // Default title
      set<title>("");
      // Default models = no selections
      set<selected,geometry>(sel::GeometryType::NO_GEOMETRY);
      set<selected,physics>(sel::PhysicsType::NO_PHYSICS);
      set<selected,position>(sel::PositionType::NO_POSITION);
      set<selected,mass>(sel::MassType::NO_MASS);
      set<selected,hydro>(sel::HydroType::NO_HYDRO);
      set<selected,energy>(sel::EnergyType::NO_ENERGY);
      set<selected,mix>(sel::MixType::NO_MIX);
      set<selected,frequency>(sel::FrequencyType::NO_FREQUENCY);
      set<selected,mixrate>(sel::MixRateType::NO_MIXRATE);
      // Default time incrementation parameters
      set<incpar,nstep>(std::numeric_limits<uint64_t>::max());
      set<incpar,term>(1.0);
      set<incpar,dt>(0.5);
      // Default number of components
      set<component,nposition>(0);
      set<component,ndensity>(0);
      set<component,nvelocity>(0);
      set<component,nscalar>(0);
      set<component,nfrequency>(0);
      set<component,npar>(1);
      // Default intervals
      set<interval,tty>(1);
      set<interval,dump>(1);
      set<interval,plot>(1);
      set<interval,pdf>(1);
      set<interval,glob>(1);
      // Default I/O parameters
      set<io,control>("");
      set<io,input>("");
      set<io,output>("out");
      set<io,pdf>("pdf");
      set<io,glob>("glob");
      set<io,stat>("stat");
      // Default beta mass model parameters
      set<param,beta,atwood>(0.5);
      // Default Dirichlet mix model parameters
      set<param,dirichlet,b>(std::vector<real>());
      set<param,dirichlet,S>(std::vector<real>());
      set<param,dirichlet,kappa>(std::vector<real>());
      // Default generalized Dirichlet mix model parameters
      set<param,gendirichlet,b>(std::vector<real>());
      set<param,gendirichlet,S>(std::vector<real>());
      set<param,gendirichlet,kappa>(std::vector<real>());
      set<param,gendirichlet,c>(std::vector<real>());
      // Default gamma mix model parameters
      set<param,gamma,c1>(0.5);
      set<param,gamma,c2>(0.73);
      set<param,gamma,c3>(5.0);
      set<param,gamma,c4>(0.25);
      // Default simplified Langevin hydro model parameters
      set<param,slm,c0>(2.1);
      // Default generalized Langevin hydro model parameters
      set<param,slm,c0>(2.1);
      // Default requested statistics
      set<stat>(std::vector<Product>());
    }

    //! Destructor
    ~InputDeck() noexcept override = default;

    //! Return total number of particle properties
    uint32_t nprop() const noexcept {
      using namespace ctr;
      return get<component,nposition>() +
             get<component,ndensity>() +
             get<component,nvelocity>() +
             get<component,nscalar>() +
             get<component,nfrequency>();
    }

    //! Return position offset
    int positionOffset() const noexcept {
      return 0;
    }
    //! Return density offset
    int densityOffset() const noexcept {
      using namespace ctr;
      return get<component,nposition>();
    }
    //! Return velocity offset
    int velocityOffset() const noexcept {
      using namespace ctr;
      return get<component,nposition>() +
             get<component,ndensity>();
    }
    //! Return scalar offset
    int scalarOffset() const noexcept {
      using namespace ctr;
      return get<component,nposition>() +
             get<component,ndensity>() +
             get<component,nvelocity>();
    }

    //! Return offset for term::quantity
    int termOffset(ctr::Quantity q) const noexcept {
      using namespace ctr;
      int offset = 0;
      if (q == Quantity::SCALAR)
        offset += get<component,nvelocity>();
      if (q == Quantity::VELOCITY_Z)
        offset += get<component,nvelocity>();
      if (q == Quantity::VELOCITY_Y)
        offset += get<component,nvelocity>();
      if (q == Quantity::VELOCITY_X)
        offset += get<component,ndensity>();
      if (q == Quantity::DENSITY)
        offset += 3 * get<component,nposition>();
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
static const InputDeck QuinoaDefaults;

} // namespace quinoa

#endif // QuinoaInputDeck_h
