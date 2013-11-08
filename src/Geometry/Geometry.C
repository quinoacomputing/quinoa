//******************************************************************************
/*!
  \file      src/Geometry/Geometry.C
  \author    J. Bakosi
  \date      Thu 07 Nov 2013 09:59:32 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry base
  \details   Geometry base
*/
//******************************************************************************

#include <Geometry.h>

using namespace quinoa;

Geometry::Geometry(const Base& base) : m_base(base)
//******************************************************************************
//  Constructor
//! \param[in]  base     Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  //! Echo information on geometry to be created
  echo();
}

void
Geometry::echo()
//******************************************************************************
//  Echo information on geometry
//! \author J. Bakosi
//******************************************************************************
{
  const QuinoaPrint& print = m_base.print;
  const ctr::InputDeck& control = m_base.control;

  print.Section<ctr::Geometry, ctr::selected, ctr::geometry>();

  print.subsection("I/O filenames");
  print.item( "Input", control.get< ctr::cmd, ctr::io, ctr::input >() );
  print.item( "Output", control.get< ctr::cmd, ctr::io, ctr::output >() );
  print.endsubsection();
}
