//******************************************************************************
/*!
  \file      src/Base/Base.h
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 04:26:27 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Collection of essentials
  \details   Collection of essentials
*/
//******************************************************************************
#ifndef Base_h
#define Base_h

#include <Print.h>
#include <Paradigm.h>
#include <Memory.h>
#include <Control.h>
#include <Timer.h>

namespace quinoa {

//! Base: collection of essentials
struct Base {
  QuinoaPrint& print;
  Paradigm& paradigm;
  Memory& memory;
  QuinoaControl& control;
  Timer& timer;

  //! Initializer constructor
  Base(QuinoaPrint& prn,
       Paradigm& par,
       Memory& mem,
       QuinoaControl& ctr,
       Timer& tim) : print(prn),
                     paradigm(par),
                     memory(mem),
                     control(ctr),
                     timer(tim) {}
};


} // namespace quinoa

#endif // Base_h
