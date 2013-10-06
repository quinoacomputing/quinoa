//******************************************************************************
/*!
  \file      src/Base/RNGTestPrint.h
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 02:30:51 PM MDT
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

//! RNGTestPrint : quinoa::Print
class RNGTestPrint : public quinoa::Print {

  public:
    //! Constructor
    explicit RNGTestPrint(const InputDeck& control) : m_ctr(control) {
      IGNORE(m_ctr);
    }

    //! Destructor
    ~RNGTestPrint() noexcept override = default;

  private:
    //! Don't permit copy constructor
    RNGTestPrint(const RNGTestPrint&) = delete;
    //! Don't permit copy assigment
    RNGTestPrint& operator=(const RNGTestPrint&) = delete;
    //! Don't permit move constructor
    RNGTestPrint(RNGTestPrint&&) = delete;
    //! Don't permit move assigment
    RNGTestPrint& operator=(RNGTestPrint&&) = delete;

    const InputDeck& m_ctr;         //!< Parsed control
};

} // rngtest::

#endif // RNGTestPrint_h
