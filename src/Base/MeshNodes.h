//******************************************************************************
/*!
  \file      src/Base/MeshNodes.h
  \author    J. Bakosi
  \date      Sun 31 Jan 2016 07:16:21 AM MST
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     MeshNodes used to store mesh node data.
  \details   MeshNodes used to store data at mesh nodes as a specialization of
    DataLayout. See also Base/DataLayout.h and the rationale discussed in the
    [design](layout.html) document.
*/
//******************************************************************************
#ifndef MeshNodes_h
#define MeshNodes_h

#include "Config.h"
#include "DataLayout.h"

namespace tk {

//! Select data layout policy for mesh node properties at compile-time
#if   defined MESHNODE_DATA_LAYOUT_AS_MESHNODE_MAJOR
using MeshNodes = DataLayout< UnkEqComp >;
#elif defined MESHNODE_DATA_LAYOUT_AS_EQUATION_MAJOR
using MeshNodes = DataLayout< EqCompUnk >;
#endif

} // tk::

#endif // MeshNodes_h
