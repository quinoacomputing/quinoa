/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/

#ifndef util_Basics_hpp
#define util_Basics_hpp

namespace phdmesh {

/**
 * \{
 * \author H. Carter Edwards  <hcedwar@sandia.gov>
 * \date   November 2006
 */

//----------------------------------------------------------------------
/** \class StaticAssert
 *  \brief  Compiler-enforced value of 'expression == true'
 *
 *  If 'expression == true' then the specialization defines
 *  - <b> enum { OK = true };                </b>
 *  - <b> static bool ok() { return true ; } </b>
 */
template< bool expression > struct StaticAssert {};

template<> struct StaticAssert<true> {
  enum { OK = true };
  static bool ok() { return true ; }
};

//----------------------------------------------------------------------
// Selection of an integer type based upon sign and size

#ifndef DOXYGEN_COMPILE

template<
  unsigned size_char  = sizeof( unsigned char  ) ,
  unsigned size_short = sizeof( unsigned short ) ,
  unsigned size_int   = sizeof( unsigned int   ) ,
  unsigned size_long  = sizeof( unsigned long  ) ,
  unsigned size_ptr   = sizeof( void * ) >
struct IntegerFundamentalTypes ;

template<>
struct IntegerFundamentalTypes<1,2,4,8,4> {
  typedef signed   char   int8_type ;
  typedef unsigned char   uint8_type ;
  typedef signed   short  int16_type ;
  typedef unsigned short  uint16_type ;
  typedef signed   int    int32_type ;
  typedef unsigned int    uint32_type ;
  typedef signed   long   int64_type ;
  typedef unsigned long   uint64_type ;
  typedef signed   int    intptr_type ;
  typedef unsigned int    uintptr_type ;
};

template<>
struct IntegerFundamentalTypes<1,2,4,8,8> {
  typedef signed   char   int8_type ;
  typedef unsigned char   uint8_type ;
  typedef signed   short  int16_type ;
  typedef unsigned short  uint16_type ;
  typedef signed   int    int32_type ;
  typedef unsigned int    uint32_type ;
  typedef signed   long   int64_type ;
  typedef unsigned long   uint64_type ;
  typedef signed   long   intptr_type ;
  typedef unsigned long   uintptr_type ;
};

template<>
struct IntegerFundamentalTypes<1,2,4,4,4> {
  typedef signed   char       int8_type ;
  typedef unsigned char       uint8_type ;
  typedef signed   short      int16_type ;
  typedef unsigned short      uint16_type ;
  typedef signed   int        int32_type ;
  typedef unsigned int        uint32_type ;
  typedef signed   long long  int64_type ;
  typedef unsigned long long  uint64_type ;
  typedef signed   long       intptr_type ;
  typedef unsigned long       uintptr_type ;
};

template<>
struct IntegerFundamentalTypes<1,2,4,4,8> {
  typedef signed   char       int8_type ;
  typedef unsigned char       uint8_type ;
  typedef signed   short      int16_type ;
  typedef unsigned short      uint16_type ;
  typedef signed   int        int32_type ;
  typedef unsigned int        uint32_type ;
  typedef signed   long long  int64_type ;
  typedef unsigned long long  uint64_type ;
  typedef signed   long long  intptr_type ;
  typedef unsigned long long  uintptr_type ;
};

#endif /* DOXYGEN_COMPILE */

typedef IntegerFundamentalTypes<>::int8_type     int8_type ;
typedef IntegerFundamentalTypes<>::uint8_type    uint8_type ;
typedef IntegerFundamentalTypes<>::int16_type    int16_type ;
typedef IntegerFundamentalTypes<>::uint16_type   uint16_type ;
typedef IntegerFundamentalTypes<>::int32_type    int32_type ;
typedef IntegerFundamentalTypes<>::uint32_type   uint32_type ;
typedef IntegerFundamentalTypes<>::int64_type    int64_type ;
typedef IntegerFundamentalTypes<>::uint64_type   uint64_type ;
typedef IntegerFundamentalTypes<>::intptr_type   intptr_type ;
typedef IntegerFundamentalTypes<>::uintptr_type  uintptr_type ;

/** \} */

} // namespace phdmesh

#endif

