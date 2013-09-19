//******************************************************************************
/*!
  \file      src/Control/QuinoaControl.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:27:04 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa control
  \details   Quinoa control
*/
//******************************************************************************
#ifndef QuinoaControl_h
#define QuinoaControl_h

#include <string>
#include <iostream>
#include <sstream>
#include <limits>

#include <QuinoaTypes.h>
#include <QuinoaControlTypes.h>
#include <Control.h>
#include <Exception.h>
#include <Option.h>

namespace quinoa {

//! QuinoaControl : Control<specialized to Quinoa>, see QuinoaControlTypes.h
class QuinoaControl :
  public Control< // tag           type
                  ctr::title,      std::string,
                  ctr::selected,   ctr::selects,
                  ctr::incpar,     ctr::incpars,
                  ctr::component,  ctr::components,
                  ctr::interval,   ctr::intervals,
                  ctr::io,         ctr::ios,
                  ctr::param,      ctr::parameters,
                  ctr::stats,      std::vector<ctr::Product> > {

  public:
    //! Constructor: set all defaults, see QuinoaControlTypes.h
    QuinoaControl() {
      using namespace ctr;
      // Default title
      set<title>("");
      // Default options
      set<selected,geometry>(select::GeometryType::NO_GEOMETRY);
      set<selected,physics>(select::PhysicsType::NO_PHYSICS);
      set<selected,position>(select::PositionType::NO_POSITION);
      set<selected,mass>(select::MassType::NO_MASS);
      set<selected,hydro>(select::HydroType::NO_HYDRO);
      set<selected,energy>(select::EnergyType::NO_ENERGY);
      set<selected,mix>(select::MixType::NO_MIX);
      set<selected,frequency>(select::FrequencyType::NO_FREQUENCY);
      set<selected,mixrate>(select::MixRateType::NO_MIXRATE);
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
      set<io,stats>("stat");
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
      set<stats>(std::vector<Product>());
    }

    //! Destructor
    ~QuinoaControl() noexcept override = default;

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
        offset += NCOMP_POS * get<component,nposition>();
      return offset;
    }

  private:
    //! Don't permit copy constructor
    QuinoaControl(const QuinoaControl&) = delete;
    //! Don't permit copy assigment
    QuinoaControl& operator=(const QuinoaControl&) = delete;
    //! Don't permit move constructor
    QuinoaControl(QuinoaControl&&) = delete;
    //! Don't permit move assigment
    QuinoaControl& operator=(QuinoaControl&&) = delete;
};

} // namespace quinoa

#endif // QuinoaControl_h
