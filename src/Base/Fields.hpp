// *****************************************************************************
/*!
  \file      src/Base/Fields.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Fields used to store data associated to mesh entities.
  \details   Fields used to store data associated to mesh entities, e.g., nodes,
    elemewnts, faces, etc., as a specialization of tk::Data. See also
    Base/Data.h and the rationale discussed in the [design](layout.html)
    document.
*/
// *****************************************************************************
#ifndef Fields_h
#define Fields_h

#include "QuinoaConfig.hpp"
#include "Data.hpp"

namespace tk {

//! Select data layout policy for mesh node properties at compile-time
#if   defined FIELD_DATA_LAYOUT_AS_FIELD_MAJOR
using Fields = Data< UnkEqComp >;
#elif defined FIELD_DATA_LAYOUT_AS_EQUATION_MAJOR
using Fields = Data< EqCompUnk >;
#endif

} // tk::

#endif // Fields_h
