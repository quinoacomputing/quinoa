//******************************************************************************
/*!
  \file      src/Mesh/GmshTxtMeshReader2D.C
  \author    J. Bakosi
  \date      Sat 06 Oct 2012 07:44:51 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh reader class definition
  \details   Gmsh mesh reader class definition
*/
//******************************************************************************

#include <sstream>
#include <iostream>

#include <GmshTxtMeshReader2D.h>

using namespace Quinoa;

void
GmshTxtMeshReader2D::read()
//******************************************************************************
//  Public interface for read 2D txt Gmsh mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Read in txt Gmsh mesh into STL containers
  GmshTxtMeshReader::read();
  // Pack STL containers into memory entries
  //pack();
}

void
GmshTxtMeshReader2D::pack()
//******************************************************************************
//  Pack mesh data from STL to memory entries
//! \author J. Bakosi
//******************************************************************************
{
//   // Count up number of elements + connectivity entries in all element sets for
//   // each element types
//   Int numLines = 0;
//   Int numTriangles = 0;
//   for (Int k=1; k<=m_elemsets; k++) {
//     stringstream ess;
//     ess << k;
//     Int* elmtype = m_memory->getPtr<Int>(ELEMTYPE_NAME+ess.str());
//     size_t num = m_memory->getNumber(ELEMID_NAME+ess.str());
//     for (Int i=0; i<num; i++) {
//       switch (elmtype[i]) {
//         case 1: numLines++; break;
//         case 2: numTriangles++; break;
//       }
//     }
//   }
// 
//   // Create new memory entry for a single packed element set storing IDs + nodes
//   // for line elements
//   m_mesh->m_connLine = m_memory->newEntry(numLines*2,
//                                           ValType::INT,
//                                           VarType::SCALAR,
//                                           LINECONN_NAME);
//   // for triangle elements
//   m_mesh->m_connTri = m_memory->newEntry(numTriangles*3,
//                                          ValType::INT,
//                                          VarType::SCALAR,
//                                          TRICONN_NAME);
// 

//    Int* elmtype = m_memory->getPtr<Int>(ELEMTYPE_NAME+ess.str());
//    pair<size_t,Int*> elmEntry = m_memory->getNumPtr<Int>(ELEMID_NAME+ess.str());

//     // elm-number elm-type number-of-tags < tag > ... node-number-list
//     for (size_t i=0; i<num; i++) {
//       cout << "  " << element[i] << " " << elmtype[i] << " {";
//       copy(m_tag[i].begin(), m_tag[i].end()-1,ostream_iterator<Int>(cout,", "));
//       cout << m_tag[i].back()-1 << "} {";
//       copy(m_elem[i].begin(),m_elem[i].end()-1,ostream_iterator<Int>(cout,", "));
//       cout << m_elem[i].back()-1 << "}" << endl;
//     }
//   }

}
