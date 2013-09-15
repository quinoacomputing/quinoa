//******************************************************************************
/*!
  \file      src/Base/Base.h
  \author    J. Bakosi
  \date      Sat 14 Sep 2013 07:35:52 PM MDT
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
  const Print& print;
  const Paradigm& paradigm;
  const Memory& memory;
  const QuinoaControl& control;
  const Timer& timer;

  //! Initializer constructor
  Base(const Print& prn,
       const Paradigm& par,
       const Memory& mem,
       const QuinoaControl& ctr,
       const Timer& tim) : print(prn),
                           paradigm(par),
                           memory(mem),
                           control(ctr),
                           timer(tim) {}
};


} // namespace quinoa

#endif // Base_h
