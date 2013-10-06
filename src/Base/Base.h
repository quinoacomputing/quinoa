//******************************************************************************
/*!
  \file      src/Base/Base.h
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 02:24:59 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Collection of essentials
  \details   Collection of essentials
*/
//******************************************************************************
#ifndef Base_h
#define Base_h

#include <QuinoaPrint.h>
#include <RNGTestPrint.h>
#include <Paradigm.h>
#include <Control.h>
#include <Timer.h>

namespace quinoa {

//! Base: collection of essentials
struct Base {
  QuinoaPrint& print;
  Paradigm& paradigm;
  InputDeck& control;
  Timer& timer;

  //! Initializer constructor
  Base(QuinoaPrint& print,
       Paradigm& paradigm,
       InputDeck& control,
       Timer& timer) : print(print),
                       paradigm(paradigm),
                       control(control),
                       timer(timer) {}
};

} // quinoa::

namespace rngtest {

//! Base: collection of essentials
struct Base {
  RNGTestPrint& print;
  quinoa::Paradigm& paradigm;
  InputDeck& control;
  quinoa::Timer& timer;

  //! Initializer constructor
  Base(RNGTestPrint& print,
       quinoa::Paradigm& paradigm,
       InputDeck& control,
       quinoa::Timer& timer) : print(print),
                               paradigm(paradigm),
                               control(control),
                               timer(timer) {}
};

} // rngtest::

#endif // Base_h
