//******************************************************************************
/*!
  \file      src/Control/QuinoaControl.h
  \author    J. Bakosi
  \date      Sun 08 Sep 2013 07:55:06 PM MDT
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
  public Control< // tag               type
                  control::title,      std::string,
                  control::selected,   control::selects,
                  control::incpar,     control::incpars,
                  control::component,  control::components,
                  control::interval,   control::intervals,
                  control::io,         control::ios,
                  control::param,      control::parameters,
                  control::stats,      std::vector<control::Product> > {

  public:
    //! Constructor: set all defaults, see QuinoaControlTypes.h
    QuinoaControl() {
      using namespace control;
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

    //! Echo vector of vector of element names
    //! Fields of vector<vector< struct{field, name, plot} >> must exist
    //! See src/Control/ControlTypes.h for the definitions of operator << for
    //! outputing Term and vector<Term>, and operator <<= for outputing
    //! requested (i.e., plotted) Term
    template<typename... tags>
    void echoVecVecNames(const std::string& msg, bool req = false) const {
      std::cout << "   - " << msg << ": {";
      if (req) {
        for (auto& v : get<tags...>()) {
          std::cout <<= v;
        }
      } else {
        for (auto& v : get<tags...>()) {
          std::cout << v;
        }
      }
      std::cout << " }" << std::endl;
    }

    //! Echo vector of Option.names
    template<class OptionType, typename... tags>
    void echoVecOptName(const std::string& msg) const {
      control::Option<OptionType> opt;
      std::cout << "   - " << msg << ": {";
      for (auto& v : get<tags...>()) {
        std::cout << " " << opt.name(v);
      }
      std::cout << " }" << std::endl;
    }

    //! Return total number of particle properties
    uint32_t nprop() const noexcept {
      using namespace control;
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
      using namespace control;
      return get<component,nposition>();
    }
    //! Return velocity offset
    int velocityOffset() const noexcept {
      using namespace control;
      return get<component,nposition>() +
             get<component,ndensity>();
    }
    //! Return scalar offset
    int scalarOffset() const noexcept {
      using namespace control;
      return get<component,nposition>() +
             get<component,ndensity>() +
             get<component,nvelocity>();
    }

    //! Return offset for term::quantity
    int termOffset(control::Quantity q) const noexcept {
      using namespace control;
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

    //! Error out on model configured at compile-time not matching that whose
    //! coefficients have been parsed
    template<class OptionType, class ModelType, typename Parsed>
    void matchModels(const ModelType configured) {
//       bool match = configured == get<options>();
//       if (!match) {
//         control::Option<OptionType> Model;
//         std::stringstream ss;
//         ss << "Compile-time-configured model (" << Model.name(configured)
//            << ") does not match that (" << Model.name(get<Parsed>())
//            << ") whose coefficients have been parsed in from the input file.";
//         Throw(ExceptType::FATAL, ss.str());
//       }
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
