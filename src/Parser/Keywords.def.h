//******************************************************************************
/*!
  \file      src/Parser/Keywords.def.h
  \author    J. Bakosi
  \date      Fri 25 Jan 2013 07:41:58 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Keywords
  \details   Includes all keywords recognized by the parser
*/
//******************************************************************************

#include <MainKeywords.def.h>
#include <PhysicsKeywords.def.h>
#include <ModelKeywords.def.h>
#include <TimeIntegrationKeywords.def.h>

namespace keyword {

  struct any : sor< main,
                    physics,
                    model,
                    time_integration > {};

} // namespace keyword
