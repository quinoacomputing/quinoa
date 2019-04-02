// *****************************************************************************
/*!
  \file      src/DiffEq/Position/PositionCoeffPolicy.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Position SDE coefficients policies
  \details   This file defines coefficients policy classes for the position
             SDE, defined in DiffEq/Position/Position.h. For general
             requirements on the position SDE coefficients policy classes see
             the header file.
*/
// *****************************************************************************

#include "PositionCoeffPolicy.h"

walker::PositionInstVel::PositionInstVel( std::array< tk::real, 9 >& )
// *****************************************************************************
//! Constructor: prescribe mean shear
// *****************************************************************************
{
}

walker::PositionConstShear::PositionConstShear( std::array< tk::real, 9 >& dU )
// *****************************************************************************
//! Constructor: prescribe mean shear
//! \param[in,out] dU Prescribed mean velocity gradient
// *****************************************************************************
{
  dU = {{ 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0,
          0.0, 0.0, 0.0 }};
}
