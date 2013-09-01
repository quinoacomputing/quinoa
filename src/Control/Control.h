//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Sat 31 Aug 2013 07:46:36 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Control base
  \details   Control base
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

#include <string>
#include <iostream>
#include <tuple>

namespace quinoa {

//! Control base templated on a tuple, its boolean version, and their fieldnames
template< typename Tuple, typename BoolTuple, typename Field >
class Control {

  protected:
    Tuple m_data;               //!< Data parsed
    BoolTuple m_booldata;       //!< Flags indicating if data was parsed

  public:
    //! Constructor: set defaults
    explicit Control(const Tuple& data) noexcept : m_data(data) {}

    //! Destructor
    virtual ~Control() noexcept = default;

    //! Set: load data
    void set(const Tuple& data, const BoolTuple& booldata) {
      m_data = move(data);
      m_booldata = move(booldata);
    }

    //! Get single element 'at' position
    template< Field at >
    constexpr const typename std::tuple_element<at, decltype(m_data)>::type&
    get() const noexcept {
      return std::get<at>(m_data);
    }

    //! Check if an element is set via boolean data
    template< Field at >
    constexpr bool set() const noexcept { return m_booldata[at]; }

    //! Echo element if set
    template< Field at >
    void echo(const std::string& msg) const {
      if (set<at>())
        std::cout << "   - " << msg << ": " << get<at>() << std::endl;
    }

    //! Echo vector of elements if set
    template< Field at >
    void echoVec(const std::string& msg) const {
      if (set<at>()) {
        std::cout << "   - " << msg << ": {";
        for (auto& v : get<at>()) std::cout << " " << v;
        std::cout << " }" << std::endl;
      }
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

} // namespace quinoa

#endif // Control_h
