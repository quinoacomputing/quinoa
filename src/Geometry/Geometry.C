//******************************************************************************
/*!
  \file      src/Geometry/Geometry.C
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 10:02:24 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry base
  \details   Geometry base
*/
//******************************************************************************

#include <Geometry.h>

using quinoa::Geometry;

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

  print.Section<ctr::Geometry, tag::selected, tag::geometry>();

  print.subsection("I/O filenames");
  print.item( "Input", control.get< tag::cmd, tag::io, tag::input >() );
  print.item( "Output", control.get< tag::cmd, tag::io, tag::output >() );
  print.endsubsection();
}
