//******************************************************************************
/*!
  \file      src/Base/Base.h
  \author    J. Bakosi
  \date      Mon Oct  7 14:29:08 2013
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
  tk::Paradigm& paradigm;
  ctr::InputDeck& control;
  tk::Timer& timer;

  //! Initializer constructor
  Base(QuinoaPrint& Print,
       tk::Paradigm& Paradigm,
       ctr::InputDeck& Control,
       tk::Timer& Timer) : print(Print),
                           paradigm(Paradigm),
                           control(Control),
                           timer(Timer) {}
};

} // quinoa::

namespace rngtest {

//! Base: collection of essentials
struct Base {
  RNGTestPrint& print;
  tk::Paradigm& paradigm;
  ctr::InputDeck& control;
  tk::Timer& timer;

  //! Initializer constructor
  Base(RNGTestPrint& Print,
       tk::Paradigm& Paradigm,
       ctr::InputDeck& Control,
       tk::Timer& Timer) : print(Print),
                           paradigm(Paradigm),
                           control(Control),
                           timer(Timer) {}
};

} // rngtest::

#endif // Base_h
