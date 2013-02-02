//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Sat 02 Feb 2013 08:32:22 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Main control category
  \details   Main control catgeory
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

#include<string>

#include<QuinoaTypes.h>
#include<Type.h>

namespace Quinoa {

//! Control base
class Control {

  private:
    // Data parsed
    control::Bundle m_data;

  public:
    //! Constructor
    Control();

    //! Destructor
    ~Control() = default;

    // Set all data in one step by deep-move of whole tuple
    void set(const control::Bundle& stack) { m_data = move(stack); }

    //! Get single element at position
    template< control::BundlePosition at >
    const typename std::tuple_element<at, decltype(m_data)>::type& get() {
      return std::get<at>(m_data);
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
