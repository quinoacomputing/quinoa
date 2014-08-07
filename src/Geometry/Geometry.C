//******************************************************************************
/*!
  \file      src/Geometry/Geometry.C
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 05:09:43 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Geometry base
  \details   Geometry base
*/
//******************************************************************************

#include <Geometry.h>

using quinoa::Geometry;

Geometry::Geometry()
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
//   const QuinoaPrint& print = m_base.print;
//   const ctr::InputDeck& control = m_base.control;
// 
//   print.Section<ctr::Geometry, tag::selected, tag::geometry>();
// 
//   print.subsection("I/O filenames");
//   print.item( "Input", control.get< tag::cmd, tag::io, tag::input >() );
//   print.item( "Output", control.get< tag::cmd, tag::io, tag::output >() );
//   print.endsubsection();
}
