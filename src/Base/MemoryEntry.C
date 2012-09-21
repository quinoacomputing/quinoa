//******************************************************************************
/*!
  \file      src/Base/MemoryEntry.C
  \author    J. Bakosi
  \date      Fri 21 Sep 2012 12:34:25 PM MDT
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
//  \return  One line info of a memory entry
//! \author  J. Bakosi
//******************************************************************************
{
  stringstream ss;
  ss << "  " << setw(EntryWidth[0]) << m_name
     << "  " << setw(EntryWidth[1]) << m_number
     << "  " << setw(EntryWidth[2]) << ValName[static_cast<Int>(m_value)]
     << "  " << setw(EntryWidth[3]) << SizeOf[static_cast<Int>(m_value)]
     << "  " << setw(EntryWidth[4]) << VarTypeName[static_cast<Int>(m_variable)]
     << "  " << setw(EntryWidth[5]) << m_bytes
     << "  " << setw(EntryWidth[6]) << (m_plot == true ? "true" : "false")
     << "  " << setw(EntryWidth[7]) << (m_restart == true ? "true" : "false")
     << "  " << setw(EntryWidth[8]) << m_ptr
     << endl;
  return ss.str();
}
