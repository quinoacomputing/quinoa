//******************************************************************************
/*!
  \file      src/Base/MemoryEntry.C
  \author    J. Bakosi
  \date      Wed 12 Sep 2012 01:58:02 AM KST
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
     << "  " << setw(EntryWidth[2]) << ValueName[m_value]
     << "  " << setw(EntryWidth[3]) << VariableTypeName[m_variable]
     << "  " << setw(EntryWidth[4]) << m_nbytes
     << "  " << setw(EntryWidth[5]) << (m_plot == true ? "true" : "false")
     << "  " << setw(EntryWidth[6]) << (m_restart == true ? "true" : "false")
     << "  " << setw(EntryWidth[7]) << m_ptr
     << endl;
  return ss.str();
}
