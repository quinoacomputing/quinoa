/*
// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include <assert.h> 

#include "RTOp_obj_null_vtbl.h"
#include "RTOp_obj_free_free.h"

/* Local function definitions */

static int get_obj_type_num_entries(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void* instance_data
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
  *num_values  = 0;
  *num_indexes = 0;
  *num_chars   = 0;
  return 0;
}

static int obj_create( const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data, void** obj )
{
  *obj = NULL;
  return 0;
}

static int obj_reinit( const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data, void* obj )
{
  return 0;
}

static int extract_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void *       instance_data
  ,void *             obj
  ,int                num_values
  ,RTOp_value_type    value_data[]
  ,int                num_indexes
  ,RTOp_index_type    index_data[]
  ,int                num_chars
  ,RTOp_char_type     char_data[]
  )
{
  assert( obj == NULL );
  assert( num_values  == 0 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
  return 0;
}

static int load_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void *            instance_data
  ,int                     num_values
  ,const RTOp_value_type   value_data[]
  ,int                     num_indexes
  ,const RTOp_index_type   index_data[]
  ,int                     num_chars
  ,const RTOp_char_type    char_data[]
  ,void **                 obj
  )
{
  assert( obj == NULL );
  assert( num_values  == 0 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
  return 0;
}

const struct RTOp_obj_type_vtbl_t   RTOp_obj_null_vtbl =
{
  get_obj_type_num_entries
  ,obj_create
  ,obj_reinit
  ,RTOp_obj_free_free
  ,extract_state
  ,load_state
};
