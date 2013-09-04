//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Tue Sep  3 07:34:48 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Control base
  \details   Control base
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

#include <string>
#include <iostream>

#include <TaggedTuple.h>

namespace quinoa {

//! Control : tagged_tuple
template< typename... Ts >
class Control : public tagged_tuple< Ts... > {

  public:
    //! Constructor: set defaults
    explicit Control() = default;

    //! Destructor
    virtual ~Control() noexcept = default;

//     //! Set: load data
//     void set(const Tuple& data) {
//       m_data = move(data);
//     }

//     //! Const-ref accessor to single element based on tag and subtag
//     template< typename tag, typename subtag >
//     constexpr const auto&
//     get() const noexcept {
//       return this->template get<tag>().template get<subtag>();
//     }

//     //! Check if an element is set
//     template< typename tag >
//     //constexpr bool set() const noexcept { return m_booldata[at]; }
//     constexpr bool set() const noexcept { return true; }

    //! Echo element if set
    template< typename tag, typename subtag >
    void echo(const std::string& msg) const {
      //if (set<tag>())
        std::cout << "   - " << msg << ": "
                  << this->template get<tag>().template get<subtag>()
                  << std::endl;
    }

//     //! Echo vector of elements if set
//     template< typename tag >
//     void echoVec(const std::string& msg) const {
//       if (set<tag>()) {
//         std::cout << "   - " << msg << ": {";
//         for (auto& v : get<tag>()) std::cout << " " << v;
//         std::cout << " }" << std::endl;
//       }
//     }

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
