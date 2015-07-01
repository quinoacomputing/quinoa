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

#ifndef util_NumericEnum_hpp
#define util_NumericEnum_hpp

#include <complex>
#include <util/Basics.hpp>
#include <util/TypeList.hpp>

namespace phdmesh {

/** \defgroup util_numeric_enum  Enumeration of Numeric Types
 *  \brief   A TypeList based enumeration of numeric types.
 *  \author  H. Carter Edwards  <hcedwar@sandia.gov>
 *  \date    October 2007
 *
 *  A TypeList is used instead of declaring an <b> enum </b>
 *  so that the integer values are clearly defined and extension of the
 *  list is robust.
 *
 *  An enumeration value for a type is given by:
 *  <b> NumericEnum< Type >::value </b>.
 *  Thus a switch statement for a numeric type is formed as follows.
 *
 *  switch( my_value ) { <br>
 *  case NumericEnum< int >::value : ... ; break ; <br>
 *  case NumericEnum< float >::value : ... ; break ; <br>
 *  case NumericEnum< double >::value : ... ; break ; <br>
 *  default: ... ; <br>
 *  } <br>
 *
 *  \{
 */

#ifndef DOXYGEN_COMPILE

typedef TypeList<          void ,
        TypeList< signed   char ,
        TypeList< unsigned char ,
        TypeList< signed   short ,
        TypeList< unsigned short ,
        TypeList< signed   int ,
        TypeList< unsigned int ,
        TypeList< signed   long ,
        TypeList< unsigned long ,
        TypeList<          float ,
        TypeList<          double ,
        TypeList<          std::complex<float> ,
        TypeList<          std::complex<double> ,

        TypeList<          void * ,
        TypeList< signed   char * ,
        TypeList< unsigned char * ,
        TypeList< signed   short * ,
        TypeList< unsigned short * ,
        TypeList< signed   int * ,
        TypeList< unsigned int * ,
        TypeList< signed   long * ,
        TypeList< unsigned long * ,
        TypeList<          float * ,
        TypeList<          double * ,
        TypeList<          std::complex<float> * ,
        TypeList<          std::complex<double> * ,

        TypeListEnd > > > > > > > > > > > > >
                    > > > > > > > > > > > > > NumericTypeList ;

#endif /* DOXYGEN_COMPILE */

template<typename Type = void> struct NumericEnum ;

/** \brief  Map the integer value associated with a numeric scalar type
 *          to a text name or byte size.
 */
template<>
struct NumericEnum<void> {

  enum { length  = TypeListLength<NumericTypeList>::value };
  enum { minimum = 1 };
  enum { maximum = length - 1 };

  /** \brief  Text name Type where where ordinal = NumericEnum<Type>::value */
  static const char * name( unsigned ordinal );

  /** \brief  sizeof(Type) where ordinal = NumericEnum<Type>::value */
  static unsigned     size( unsigned ordinal );

  enum { value = TypeListIndex< NumericTypeList , void>::value };

private:

  enum { OK = StaticAssert< TypeListUnique<NumericTypeList>::value >::OK };
};

/** \brief  Map a numeric scalar Type to an integer value */
template<typename Type>
struct NumericEnum {

  /** \brief  Unique integer value for numeric scalar Type, if known.
   *          Is -1 if the type is unknown.
   */
  enum { value
#ifndef DOXYGEN_COMPILE
               = TypeListIndex<NumericTypeList,Type>::value
#endif
  };

private:
  enum { OK = StaticAssert<
         (0 < (int) value) &&
         ((int) value < (int) NumericEnum<void>::length) >::OK };
};

/** \brief  Inverse map of a numeric scalar type to an integer value */
template<unsigned Ordinal>
struct NumericType {
#ifndef DOXYGEN_COMPILE
  typedef typename TypeListAt< NumericTypeList , Ordinal >::type type ;
#endif
private:
  enum { OK = StaticAssert< Ordinal < NumericEnum<>::length >::OK };
};

/** \} */

} // namespace phdmesh

#endif

