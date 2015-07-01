/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2005 Lawrence Livermore National Laboratory.  Under 
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


/** \file ArrayIterator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ArrayIterator.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

ArrayIterator::~ArrayIterator()
  {}
  
void ArrayIterator::restart()
  { mIter = mArray.begin(); }

Mesh::EntityHandle ArrayIterator::operator*() const
  { return *mIter; }

void ArrayIterator::operator++()
  { ++mIter; }

bool ArrayIterator::is_at_end()
  { return mIter == mArray.end(); }

VertexArrayIterator::VertexArrayIterator( Mesh* meshInterface, MsqError& err )
{
  meshInterface->get_all_vertices( mArray, err );
  if (MSQ_CHKERR(err))
    mArray.clear();
  restart();
}

ElementArrayIterator::ElementArrayIterator( Mesh* meshInterface, MsqError& err )
{
  meshInterface->get_all_elements( mArray, err );
  if (MSQ_CHKERR(err))
    mArray.clear();
  restart();
}

} // namespace Mesquite
