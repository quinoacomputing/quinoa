//******************************************************************************
/*!
  \file      src/Base/Base.h
  \author    J. Bakosi
  \date      Thu Oct  3 15:43:46 2013
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
  Base(QuinoaPrint& prn,
       Paradigm& par,
       InputDeck& ctr,
       Timer& tim) : print(prn),
                     paradigm(par),
                     control(ctr),
                     timer(tim) {}
};


} // namespace quinoa

#endif // Base_h
