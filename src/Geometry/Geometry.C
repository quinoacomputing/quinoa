//******************************************************************************
/*!
  \file      src/Geometry/Geometry.C
  \author    J. Bakosi
  \date      Fri Sep 27 15:27:58 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry base
  \details   Geometry base
*/
//******************************************************************************

#include <Geometry.h>
#include <QuinoaControl.h>

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
  m_base.print.section<sel::Geometry, ctr::selected, ctr::geometry>();

  m_base.print.subsection("I/O filenames");
  m_base.print.item("Input", m_base.control.get<ctr::io,ctr::input>());
  m_base.print.item("Output", m_base.control.get<ctr::io,ctr::geomoutput>());
  m_base.print.endsubsection();
}
