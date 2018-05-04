// *****************************************************************************
/*!
  \file      src/Main/QuinoaConfig.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Quinoa configuration query functions based on cmake input
  \details   Quinoa configuration query functions based on cmake input.
*/
// *****************************************************************************

#include "QuinoaConfig.h"

namespace tk {

std::string commit()
// *****************************************************************************
//! Query git commit sha1
//! \return git commit sha1
// *****************************************************************************
{
  return GIT_COMMIT;
}

} // tk::
