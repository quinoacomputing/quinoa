//******************************************************************************
/*!
  \file      src/Base/MemoryEntry.C
  \author    J. Bakosi
  \date      Mon 10 Sep 2012 06:04:31 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MemoryEntry base class definition
  \details   Memoryentry base class definition
*/
//******************************************************************************

#include <sstream>
#include <iomanip>

using namespace std;

#include <MemoryEntry.h>

using namespace Quinoa;

string
MemoryEntry::line()
//******************************************************************************
//  One-liner returning all fields
//! \author  J. Bakosi
//******************************************************************************
{
  stringstream ss;
  ss << "  " << setw(10) << m_name
     << "  " << setw(10) << m_number
     << "  " << setw(10) << m_value
     << "  " << setw(10) << m_variable
     << "  " << setw(10) << m_nbytes
     << "  " << setw(10) << m_plot
     << "  " << setw(10) << m_restart
     << "  " << setw(10) << m_ptr
     << endl;
  return ss.str();
}
