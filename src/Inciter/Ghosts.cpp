// *****************************************************************************
/*!
  \file      src/Inciter/Ghosts.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Definitions file for generating ghost data structures
  \details   Definitions file for asynchronous distributed
             ghost data structures using Charm++.
*/
// *****************************************************************************

#include "Ghosts.hpp"
#include "DerivedData.hpp"
#include "Reorder.hpp"

using inciter::Ghosts;

Ghosts::Ghosts( const CProxy_Discretization& disc,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
//  m_fd( m_inpoel, bface, tk::remap(triinpoel,Disc()->Lid()) )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
}

#include "NoWarning/ghosts.def.h"
