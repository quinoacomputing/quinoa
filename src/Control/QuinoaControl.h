//******************************************************************************
/*!
  \file      src/Control/QuinoaControl.h
  \author    J. Bakosi
  \date      Wed Sep  4 12:30:36 2013
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

#include <QuinoaTypes.h>
#include <QuinoaControlTypes.h>
#include <Control.h>
#include <Exception.h>
#include <Option.h>

namespace quinoa {

//! QuinoaControl : Control<specialized to Quinoa>, see QuinoaControlTypes.h
class QuinoaControl :
  public Control< // tag              type
                  control::title,     std::string,
                  control::selected,  control::selects,
                  control::incpar,    control::incpars,
                  control::component, control::components,
                  control::interval,  control::intervals,
                  control::io,        control::ios,
                  control::parameter, control::parameters,
                  control::statistic, control::statistics > {

  public:
    //! Constructor: set defaults
    explicit QuinoaControl() = default;

//   //! Default constructor with defaults
//   explicit options() : geometry(select::GeometryType::NO_GEOMETRY),
//                        physics(select::PhysicsType::NO_PHYSICS),
//                        position(select::PositionType::NO_POSITION),
//                        mass(select::MassType::NO_MASS),
//                        hydro(select::HydroType::NO_HYDRO),
//                        energy(select::EnergyType::NO_ENERGY),
//                        mix(select::MixType::NO_MIX),
//                        frequency(select::FrequencyType::NO_FREQUENCY),
//                        mixrate(select::MixRateType::NO_MIXRATE) {}
// };
//   //! Default constructor with defaults
//   explicit incpars() : nstep(std::numeric_limits<uint64_t>::max()),
//                        term(1.0),
//                        dt(0.5) {}
//   //! Default constructor with defaults
//   explicit components() : nposition(0),
//                           ndensity(0),
//                           nvelocity(0),
//                           nscalar(0),
//                           nfrequency(0),
//                           npar(1) {}
// 
//   //! Return sum of all components
//   uint32_t sum() const noexcept {
//     return nposition + ndensity + nvelocity + nscalar + nfrequency;
//   }
//     //! Default constructor with defaults
//   explicit intervals() : tty(1),
//                          dump(1),
//                          plot(1),
//                          pdf(1),
//                          glob(1) {}
//   //! Default constructor with defaults
//   explicit ios() : input(),
//                    output(),
//                    pdf("jpdf"),
//                    glob("glob"),
//                    stat("stat") {}

    //! Destructor
    ~QuinoaControl() noexcept override = default;

    //! Echo vector of vector of element names if set
    //! Fields of vector<vector< struct{field, name, plot} >> must exist
    //! See src/Control/ControlTypes.h for the definitions of operator << for
    //! outputing Term and vector<Term>, and operator <<= for outputing
    //! requested (i.e., plotted) Term
    template< typename tag >
    void echoVecVecNames(const std::string& msg, bool req = false) const {
      //if (set<at>()) {
        std::cout << "   - " << msg << ": {";
        if (req) for (auto& v : this->template get<tag>()) std::cout <<= v;
        else for (auto& v : this->template get<tag>()) std::cout << v;
        std::cout << " }" << std::endl;
      //}
    }

    //! Echo vector of Option.names if set
    template< typename tag, class OptionType >
    void echoVecOptName(const std::string& msg) const {
      //if (set<at>()) {
        control::Option<OptionType> opt;
        std::cout << "   - " << msg << ": {";
        for (auto& v : this->template get<tag>()) {
          std::cout << " " << opt.name(v);
        }
        std::cout << " }" << std::endl;
      //}
    }

    //! Return total number of particle properties
    uint32_t nprop() const noexcept {
      return 1;//this->template get<control::component>().sum();
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
