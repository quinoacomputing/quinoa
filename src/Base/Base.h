//******************************************************************************
/*!
  \file      src/Base/Base.h
  \author    J. Bakosi
  \date      Sat 05 Oct 2013 05:03:20 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Collection of essentials
  \details   Collection of essentials
*/
//******************************************************************************
#ifndef Base_h
#define Base_h

#include <Print.h>
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


} // namespace quinoa

#endif // Base_h
