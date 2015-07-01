/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_LAPLACIAN_COMMON_CPP
#define MSQ_LAPLACIAN_COMMON_CPP

#include "LaplacianCommon.hpp"

namespace MESQUITE_NS {

LaplacianCommon::LaplacianCommon( ObjectiveFunction* OF ) 
  : VertexMover( OF, true ),
    patchSet( 1, true ) 
  {}

LaplacianCommon::~LaplacianCommon() {}

PatchSet* LaplacianCommon::get_patch_set() { return &patchSet; }

} // namespace Mesquite

#endif
