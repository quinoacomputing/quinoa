// *****************************************************************************
/*!
  \file      src/Inciter/ConjugateGradients.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare array for distributed conjugate gradients.
  \details   Charm++ chare array for asynchronous distributed
    conjugate gradients linear solver.
*/
// *****************************************************************************

#include "ConjugateGradients.hpp"

using tk::ConjugateGradients;

ConjugateGradients::ConjugateGradients(
  std::size_t size,
  std::size_t dof,
  const std::pair< std::vector< std::size_t >,
                   std::vector< std::size_t > >& psup ) :
  m_A( dof, psup ),
  m_x( size, 0.0 ),
  m_r( size, 0.0 ),
  m_p( size, 0.0 ),
  m_q( size, 0.0 )
// *****************************************************************************
//  Constructor
//! \param[in] size Number of unknowns (rows) on this chare
// *****************************************************************************
{
}

#include "NoWarning/conjugategradients.def.h"
