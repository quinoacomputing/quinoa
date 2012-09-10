//******************************************************************************
/*!
  \file      src/Base/MemoryEntry.C
  \author    J. Bakosi
  \date      Mon 10 Sep 2012 01:54:38 PM KST
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
  ss << "  " << setw(MemoryEntryWidth[0]) << m_name
     << "  " << setw(MemoryEntryWidth[1]) << m_number
     << "  " << setw(MemoryEntryWidth[2]) << ValueName[m_value]
     << "  " << setw(MemoryEntryWidth[3]) << VariableTypeName[m_variable]
     << "  " << setw(MemoryEntryWidth[4]) << m_nbytes
     << "  " << setw(MemoryEntryWidth[5]) << (m_plot == true ? "true" : "false")
     << "  " << setw(MemoryEntryWidth[6]) << (m_restart == true ? "true" : "false")
     << "  " << setw(MemoryEntryWidth[7]) << m_ptr
     << endl;
  return ss.str();
}
