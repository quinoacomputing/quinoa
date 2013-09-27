//******************************************************************************
/*!
  \file      src/Base/Base.h
  \author    J. Bakosi
  \date      Thu 26 Sep 2013 11:11:53 PM MDT
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
  QuinoaControl& defctr;
  Timer& timer;

  //! Initializer constructor
  Base(QuinoaPrint& prn,
       Paradigm& par,
       Memory& mem,
       QuinoaControl& ctr,
       QuinoaControl& def,
       Timer& tim) : print(prn),
                     paradigm(par),
                     memory(mem),
                     control(ctr),
                     defctr(def),
                     timer(tim) {}
};


} // namespace quinoa

#endif // Base_h
