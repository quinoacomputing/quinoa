//******************************************************************************
/*!
  \file      src/Mesh/Mesh.C
  \author    J. Bakosi
  \date      Wed 12 Sep 2012 08:30:02 PM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh base class definition
  \details   Mesh base class definition
*/
//******************************************************************************

#include <iterator>
#include <iostream>
#include <sstream>

#include <Mesh.h>
#include <Memory.h>
#include <MemoryException.h>
#include <MeshException.h>

using namespace Quinoa;

Mesh::~Mesh()
//******************************************************************************
//  Destructor: free sets
//! \author J. Bakosi
//******************************************************************************
{
  // Free node sets
  if (m_entry.size()) {
    // Free all mesh sets
    MeshSet::const_iterator it;
    for (it=m_entry.begin(); it!=m_entry.end(); it++) {
      m_memory->freeEntry(*it);
    }
    // Clear container
    m_entry.clear();
  }
}

void
Mesh::reserveElem(vector< vector<Int> >::size_type n)
//******************************************************************************
//  Add new element
//! \param[in]  n  Desired new capacity to store n elements with their tags
//! \author J. Bakosi
//******************************************************************************
{
  try {
    m_elem.reserve(n);
    m_tag.reserve(n);
  } catch (bad_alloc& ba) { throw MemoryException(FATAL, BAD_ALLOC); }  
}

void
Mesh::addElem(vector<int>& nodes)
//******************************************************************************
//  Add new element
//! \param[in]  nodes  Vector of node ids (i.e. connectivity) of the new element
//! \author J. Bakosi
//******************************************************************************
{
  try {
    m_elem.push_back(nodes);
  } catch (bad_alloc& ba) { throw MemoryException(FATAL, BAD_ALLOC); }
}

void
Mesh::addElemTags(vector<Int>& tags)
//******************************************************************************
//  Add new element tags
//! \param[in]  tags  Vector of tags to be added
//! \author J. Bakosi
//******************************************************************************
{
  try {
    m_tag.push_back(tags);
  } catch (bad_alloc& ba) { throw MemoryException(FATAL, BAD_ALLOC); }
}

void
Mesh::echoElemSets()
//******************************************************************************
//  Echo element tags and connectivity in all element sets
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if there are no element sets
  if (!m_elemsets) throw MeshException(WARNING, EMPTY_SET);

  // Echo all element sets
  for (Int k=1; k<=m_elemsets; k++) {
    cout << "* Element set: " << k << endl << endl;
    stringstream ess;
    ess << k;
    Int* elmtype = m_memory->getPtr<Int>(ELEMTYPE_NAME+ess.str());
    pair<size_t,Int*> elmEntry = m_memory->getNumPtr<Int>(ELEMID_NAME+ess.str());
    size_t num = elmEntry.first;
    Int* element = elmEntry.second;

    // elm-number elm-type number-of-tags < tag > ... node-number-list
    for (size_t i=0; i<num; i++) {
      cout << "  " << element[i] << " " << elmtype[i] << " {";
      copy(m_tag[i].begin(),m_tag[i].end()-1,ostream_iterator<Int>(cout,", "));
      cout << m_tag[i].back()-1 << "} {";
      copy(m_elem[i].begin(),m_elem[i].end()-1,ostream_iterator<Int>(cout,", "));
      cout << m_elem[i].back()-1 << "}" << endl;
    }
  }
}
