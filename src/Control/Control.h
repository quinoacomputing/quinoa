//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Tue 19 Feb 2013 07:43:36 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Main control category
  \details   Main control catgeory
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

#include <string>

#include <QuinoaTypes.h>
#include <ControlTypes.h>
#include <BackAssociate.h>
#include <Defaults.h>

namespace Quinoa {

using namespace control;

//! Control base
class Control {

  private:
    Bundle m_data;              //! Data parsed
    BoolBundle m_booldata;      //! Flags indicating if data was parsed
    string m_jpdf_filename_base;//! This will be a bundle from cmd line parse

  public:
    //! Constructor
    Control();

    //! Destructor
    ~Control() = default;

    //! Set all data in one step by deep-move of whole bundle
    void set(const Bundle& stack) { m_data = move(stack); }

    //! Set all flags in one step by deep-move of whole bool bundle
    void set(const BoolBundle& boolstack) { m_booldata = move(boolstack); }

    //! Set jpdf filename base
    void set(const string& jpdf_filename_base) {
      m_jpdf_filename_base = move(jpdf_filename_base);
    }

    //! Get jpdf filename base
    const string& get_jpdf_filename_base() const { return m_jpdf_filename_base; }

    //! Get single element at position
    template< BundlePosition at >
    const typename std::tuple_element<at, decltype(m_data)>::type& get() const {
      return std::get<at>(m_data);
    }

    //! Check if an element is set during parse
    template< BundlePosition at >
    bool set() const { return m_booldata[at]; }

    //! Get physics keyword
    const std::string& physicsKeyword() const {
      return associate::PhysicsKeyword[ std::get<PHYSICS>(m_data) ];
    }
    //! Get physics name
    const std::string& physicsName() const {
      return associate::PhysicsName[ std::get<PHYSICS>(m_data) ];
    }

    //! Get hydrodynamics model keyword
    const std::string& hydroKeyword() const {
      return associate::HydroKeyword[ std::get<HYDRO>(m_data) ];
    }
    //! Get hydrodynamics model name
    const std::string& hydroName() const {
      return associate::HydroName[ std::get<HYDRO>(m_data) ];
    }

    //! Get material mix model keyword
    const std::string& mixKeyword() const {
      return associate::MixKeyword[ std::get<MIX>(m_data) ];
    }
    //! Get material mix model name
    const std::string& mixName() const {
      return associate::MixName[ std::get<MIX>(m_data) ];
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
