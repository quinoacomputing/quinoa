// *****************************************************************************
/*!
  \file      src/Inciter/AMR/Error.C
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Class for computing error estimates for mesh refinement
  \details   Class for computing error estimates for mesh refinement.
*/
// *****************************************************************************

#include "Error.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using AMR::Error;

Error::Error()
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
}

void
Error::scalar( const tk::Fields& u )
// *****************************************************************************
//  Estimate error for scalar quantity
//! \param[in] u Current solution vector
// *****************************************************************************
{
}
