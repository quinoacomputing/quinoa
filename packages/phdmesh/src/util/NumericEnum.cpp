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
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   October 2007
 */

#include <util/NumericEnum.hpp>

namespace phdmesh {

enum { OkNumericTypeList = StaticAssert<
  NumericEnum< void >                 ::value  ==  0  &&
  NumericEnum< signed   char >        ::value  ==  1  &&
  NumericEnum< unsigned char >        ::value  ==  2  &&
  NumericEnum< signed   short >       ::value  ==  3  &&
  NumericEnum< unsigned short >       ::value  ==  4  &&
  NumericEnum< signed   int >         ::value  ==  5  &&
  NumericEnum< unsigned int >         ::value  ==  6  &&
  NumericEnum< signed   long >        ::value  ==  7  &&
  NumericEnum< unsigned long >        ::value  ==  8  &&
  NumericEnum< float >                ::value  ==  9  &&
  NumericEnum< double >               ::value  == 10  &&
  NumericEnum< std::complex<float> >  ::value  == 11  &&
  NumericEnum< std::complex<double> > ::value  == 12  &&

  NumericEnum< void * >                 ::value  == 13  &&
  NumericEnum< signed   char * >        ::value  == 14  &&
  NumericEnum< unsigned char * >        ::value  == 15  &&
  NumericEnum< signed   short * >       ::value  == 16  &&
  NumericEnum< unsigned short * >       ::value  == 17  &&
  NumericEnum< signed   int * >         ::value  == 18  &&
  NumericEnum< unsigned int * >         ::value  == 19  &&
  NumericEnum< signed   long * >        ::value  == 20  &&
  NumericEnum< unsigned long * >        ::value  == 21  &&
  NumericEnum< float * >                ::value  == 22  &&
  NumericEnum< double * >               ::value  == 23  &&
  NumericEnum< std::complex<float> * >  ::value  == 24  &&
  NumericEnum< std::complex<double> * > ::value  == 25  &&

  NumericEnum<>                         ::length == 26 >::OK };

const char * NumericEnum<void>::name( unsigned t )
{
  static const char name_void[]           = "void" ;
  static const char name_schar[]          = "signed char" ;
  static const char name_uchar[]          = "unsigned char" ;
  static const char name_short[]          = "signed short" ;
  static const char name_ushort[]         = "unsigned short" ;
  static const char name_int[]            = "signed int" ;
  static const char name_uint[]           = "unsigned int" ;
  static const char name_long[]           = "signed long" ;
  static const char name_ulong[]          = "unsigned long" ;
  static const char name_float[]          = "float" ;
  static const char name_double[]         = "double" ;
  static const char name_complex_float[]  = "complex<float>" ;
  static const char name_complex_double[] = "complex<double>" ;

  static const char name_void_p[]           = "void *" ;
  static const char name_schar_p[]          = "signed char *" ;
  static const char name_uchar_p[]          = "unsigned char *" ;
  static const char name_short_p[]          = "signed short *" ;
  static const char name_ushort_p[]         = "unsigned short *" ;
  static const char name_int_p[]            = "signed int *" ;
  static const char name_uint_p[]           = "unsigned int *" ;
  static const char name_long_p[]           = "signed long *" ;
  static const char name_ulong_p[]          = "unsigned long *" ;
  static const char name_float_p[]          = "float *" ;
  static const char name_double_p[]         = "double *" ;
  static const char name_complex_float_p[]  = "complex<float> *" ;
  static const char name_complex_double_p[] = "complex<double> *" ;

  static const char name_error[] = "ERROR" ;

  static const char * name_list[] = {
    name_void ,
    name_schar , name_uchar ,
    name_short , name_ushort ,
    name_int ,   name_uint ,
    name_long ,  name_ulong ,
    name_float , name_double ,
    name_complex_float , name_complex_double ,

    name_void_p ,
    name_schar_p , name_uchar_p ,
    name_short_p , name_ushort_p ,
    name_int_p ,   name_uint_p ,
    name_long_p ,  name_ulong_p ,
    name_float_p , name_double_p ,
    name_complex_float_p , name_complex_double_p ,
    name_error };

  if ( length < t ) { t = length ; }

  return name_list[ t ];
}

unsigned NumericEnum<void>::size( unsigned t )
{
  static unsigned size_list[] = {
    0 ,
    sizeof( NumericType<  1 >::type ) ,
    sizeof( NumericType<  2 >::type ) ,
    sizeof( NumericType<  3 >::type ) ,
    sizeof( NumericType<  4 >::type ) ,
    sizeof( NumericType<  5 >::type ) ,
    sizeof( NumericType<  6 >::type ) ,
    sizeof( NumericType<  7 >::type ) ,
    sizeof( NumericType<  8 >::type ) ,
    sizeof( NumericType<  9 >::type ) ,
    sizeof( NumericType< 10 >::type ) ,
    sizeof( NumericType< 11 >::type ) ,
    sizeof( NumericType< 12 >::type ) ,

    sizeof( NumericType< 13 >::type ) ,
    sizeof( NumericType< 14 >::type ) ,
    sizeof( NumericType< 15 >::type ) ,
    sizeof( NumericType< 16 >::type ) ,
    sizeof( NumericType< 17 >::type ) ,
    sizeof( NumericType< 18 >::type ) ,
    sizeof( NumericType< 19 >::type ) ,
    sizeof( NumericType< 20 >::type ) ,
    sizeof( NumericType< 21 >::type ) ,
    sizeof( NumericType< 22 >::type ) ,
    sizeof( NumericType< 23 >::type ) ,
    sizeof( NumericType< 24 >::type ) ,
    sizeof( NumericType< 25 >::type ) };

  if ( length < t ) { t = 0 ; }

  return size_list[ t ];
}

}


