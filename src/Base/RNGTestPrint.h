//******************************************************************************
/*!
  \file      src/Base/RNGTestPrint.h
  \author    J. Bakosi
  \date      Thu 31 Oct 2013 09:54:03 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest's printer
  \details   RNGTest's printer
*/
//******************************************************************************
#ifndef RNGTestPrint_h
#define RNGTestPrint_h

#include <Macro.h>
#include <Print.h>
#include <RNGTest/InputDeck/InputDeck.h>

namespace rngtest {

//! RNGTestPrint : Print
class RNGTestPrint : public tk::Print {

  public:
    //! Constructor
    explicit RNGTestPrint(const std::unique_ptr< ctr::InputDeck >& control) :
      m_ctr(*control) {}

    //! Destructor
    ~RNGTestPrint() override = default;

//     //! Print control option: 'group : option' only if differs from its default
//     template<typename OptionType, typename... tags>
//     void item() const {
//       if (m_ctr.get<tags...>() != ctr::InputDeckDefaults.get<tags...>()) {
//         tk::Option<OptionType> opt;
//         std::cout << m_item_name_value_fmt % m_item_indent
//                                            % opt.group()
//                                            % opt.name(m_ctr.get<tags...>()[0]);
//       }
//     }

  private:
    //! Don't permit copy constructor
    RNGTestPrint(const RNGTestPrint&) = delete;
    //! Don't permit copy assigment
    RNGTestPrint& operator=(const RNGTestPrint&) = delete;
    //! Don't permit move constructor
    RNGTestPrint(RNGTestPrint&&) = delete;
    //! Don't permit move assigment
    RNGTestPrint& operator=(RNGTestPrint&&) = delete;

    const ctr::InputDeck& m_ctr;         //!< Parsed control
};

} // rngtest::

#endif // RNGTestPrint_h
