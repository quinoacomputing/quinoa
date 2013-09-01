//******************************************************************************
/*!
  \file      src/Control/QuinoaControl.h
  \author    J. Bakosi
  \date      Sun 01 Sep 2013 02:28:53 PM MDT
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
#include <Control.h>
#include <QuinoaControlTypes.h>
#include <Exception.h>
#include <Option.h>

namespace quinoa {

//! QuinoaControl : Control< specialized to Quinoa's control >
class QuinoaControl : public Control< control::Bundle,
                                      control::BoolBundle,
                                      control::BundlePosition > {

  public:
    //! Constructor: set defaults
    explicit QuinoaControl() noexcept : Control(control::defaults) {}

    //! Destructor
    ~QuinoaControl() noexcept override = default;

    //! Echo vector of vector of element names if set
    //! Fields of vector<vector< struct{field, name, plot} >> must exist
    //! See src/Control/ControlTypes.h for the definitions of operator << for
    //! outputing Term and vector<Term>, and operator <<= for outputing
    //! requested (i.e., plotted) Term
    template< control::BundlePosition at >
    void echoVecVecNames(const std::string& msg, bool req = false) const {
      if (set<at>()) {
        std::cout << "   - " << msg << ": {";
        if (req) for (auto& v : get<at>()) std::cout <<= v;
        else for (auto& v : get<at>()) std::cout << v;
        std::cout << " }" << std::endl;
      }
    }

    //! Echo vector of Option.names if set
    template< control::BundlePosition at, class OptionType >
    void echoVecOptName(const std::string& msg) const {
      if (set<at>()) {
        control::Option<OptionType> opt;
        std::cout << "   - " << msg << ": {";
        for (auto& v : get<at>()) std::cout << " " << opt.name(v);
        std::cout << " }" << std::endl;
      }
    }

    //! Return total number of particle properties
    int nprop() const noexcept {
      return std::get<control::NPOSITION>(m_data) +
             std::get<control::NDENSITY>(m_data) +
             std::get<control::NVELOCITY>(m_data) +
             std::get<control::NSCALAR>(m_data);
    }

    //! Return position offset
    int positionOffset() const noexcept {
      return 0;
    }
    //! Return density offset
    int densityOffset() const noexcept {
      return std::get<control::NPOSITION>(m_data);
    }
    //! Return velocity offset
    int velocityOffset() const noexcept {
      return std::get<control::NPOSITION>(m_data) +
             std::get<control::NDENSITY>(m_data);
    }
    //! Return scalar offset
    int scalarOffset() const noexcept {
      return std::get<control::NPOSITION>(m_data) +
             std::get<control::NDENSITY>(m_data) +
             std::get<control::NVELOCITY>(m_data);
    }

    //! Return offset for term::quantity
    int termOffset(control::Quantity q) const noexcept {
      using namespace control;
      int offset = 0;
      if (q == Quantity::SCALAR)     offset += std::get<NVELOCITY>(m_data);
      if (q == Quantity::VELOCITY_Z) offset += std::get<NVELOCITY>(m_data);
      if (q == Quantity::VELOCITY_Y) offset += std::get<NVELOCITY>(m_data);
      if (q == Quantity::VELOCITY_X) offset += std::get<NDENSITY>(m_data);
      if (q == Quantity::DENSITY)
        offset += NCOMP_POS * std::get<NPOSITION>(m_data);
      return offset;
    }

    //! Error out on model configured at compile-time not matching that whose
    //! coefficients have been parsed
    template<class OptionType, class ModelType, control::BundlePosition Parsed>
    void matchModels(const ModelType configured) {
      bool match = configured == get<Parsed>();
      if (!match) {
        control::Option<OptionType> Model;
        std::stringstream ss;
        ss << "Compile-time-configured model (" << Model.name(configured)
           << ") does not match that (" << Model.name(get<Parsed>())
           << ") whose coefficients have been parsed in from the input file.";
        Throw(ExceptType::FATAL, ss.str());
      }
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
