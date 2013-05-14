//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Mon 13 May 2013 08:45:27 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Main control category
  \details   Main control catgeory
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

#include <string>
#include <iostream>

#include <QuinoaTypes.h>
#include <ControlTypes.h>
#include <BackAssociate.h>
#include <Defaults.h>

namespace Quinoa {

//! Control base
class Control {

  private:
    control::Bundle m_data;        //! Data parsed
    control::BoolBundle m_booldata;//! Flags indicating if data was parsed

  public:
    //! Constructor
    explicit Control() noexcept : m_data(control::DEFAULTS) {}

    //! Destructor
    ~Control() noexcept = default;

    //! Set all data in one step by deep-move of whole bundle
    void set(const control::Bundle& stack) { m_data = move(stack); }

    //! Set all flags in one step by deep-move of whole bool bundle
    void set(const control::BoolBundle& boolstack) {
      m_booldata = move(boolstack);
    }

    //! Get single element at position
    template< control::BundlePosition at >
    const typename std::tuple_element<at, decltype(m_data)>::type& get()
    const noexcept {
      return std::get<at>(m_data);
    }

    //! Check if an element is set during parse
    template< control::BundlePosition at >
    bool set() const noexcept { return m_booldata[at]; }

    //! Echo element if set
    template< control::BundlePosition at >
    void echo(const std::string& msg) const {
      if (set<at>()) cout << "   - " << msg << ": " << get<at>() << endl;
    }

    //! Echo vector of elements if set
    template< control::BundlePosition at >
    void echoVec(const std::string& msg) const {
      if (set<at>()) {
        cout << "   - " << msg << ": {";
        for (auto& v : get<at>()) cout << " " << v;
        cout << " }" << endl;
      }
    }

    //! Echo vector of vector of element names if set.
    //! Fields of vector<vector< struct{field,name,plot} >> must exist.
    //! See src/Control/ControlTypes.h for the definitions of operator << for
    //! outputing Term and vector<Term>, and operator <<= for outputing
    //! requested (i.e., plotted) Term.
    template< control::BundlePosition at >
    void echoVecVecNames(const std::string& msg, bool req = false) const {
      if (set<at>()) {
        cout << "   - " << msg << ": {";
        if (req) for (auto& v : get<at>()) cout <<= v;
        else for (auto& v : get<at>()) cout << v;
        cout << " }" << endl;
      }
    }

    //! Get physics keyword
    const std::string& physicsKeyword() const noexcept {
      return associate::PhysicsKeyword[ std::get<control::PHYSICS>(m_data) ];
    }
    //! Get physics name
    const std::string& physicsName() const noexcept {
      return associate::PhysicsName[ std::get<control::PHYSICS>(m_data) ];
    }

    //! Get mass model keyword
    const std::string& massKeyword() const noexcept {
      return associate::MassKeyword[ std::get<control::MASS>(m_data) ];
    }
    //! Get mass model name
    const std::string& massName() const noexcept {
      return associate::MassName[ std::get<control::MASS>(m_data) ];
    }

    //! Get hydrodynamics model keyword
    const std::string& hydroKeyword() const noexcept {
      return associate::HydroKeyword[ std::get<control::HYDRO>(m_data) ];
    }
    //! Get hydrodynamics model name
    const std::string& hydroName() const noexcept {
      return associate::HydroName[ std::get<control::HYDRO>(m_data) ];
    }

    //! Get material mix model keyword
    const std::string& mixKeyword() const noexcept {
      return associate::MixKeyword[ std::get<control::MIX>(m_data) ];
    }
    //! Get material mix model name
    const std::string& mixName() const noexcept {
      return associate::MixName[ std::get<control::MIX>(m_data) ];
    }

  private:
    //! Don't permit copy constructor
    Control(const Control&) = delete;
    //! Don't permit copy assigment
    Control& operator=(const Control&) = delete;
    //! Don't permit move constructor
    Control(Control&&) = delete;
    //! Don't permit move assigment
    Control& operator=(Control&&) = delete;
};

} // namespace Quinoa

#endif // Control_h
