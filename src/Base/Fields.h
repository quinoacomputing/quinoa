// *****************************************************************************
/*!
  \file      src/Base/Fields.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Fields used to store data associated to mesh entities.
  \details   Fields used to store data associated to mesh entities, e.g., nodes,
    elemewnts, faces, etc., as a specialization of tk::Data. See also
    Base/Data.h and the rationale discussed in the [design](layout.html)
    document.
*/
// *****************************************************************************
#ifndef Fields_h
#define Fields_h

#include "QuinoaConfig.h"
#include "Data.h"

namespace tk {

//! Select data layout policy for mesh node properties at compile-time
#if   defined FIELD_DATA_LAYOUT_AS_FIELD_MAJOR
using Fields = Data< UnkEqComp >;
#elif defined FIELD_DATA_LAYOUT_AS_EQUATION_MAJOR
using Fields = Data< EqCompUnk >;
#endif

} // tk::

#endif // Fields_h
