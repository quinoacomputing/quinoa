#!/usr/bin/python

c_head = """
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>
#include "H5Part.h"
#include "H5PartErrors.h"
#include "H5PartPrivate.h"

#include "H5BlockTypes.h"
#include "H5Block.h"
#include "H5BlockPrivate.h"
#include "H5BlockErrors.h"
"""

h_head = """
#ifndef _H5BLOCK_READWRITE_H_
#define _H5BLOCK_READWRITE_H_

#ifdef __cplusplus
extern "C" {
#endif

"""

h_tail = """

#ifdef __cplusplus
}
#endif

#endif
"""

fc_head = """
#include "H5Part.h"
#include "H5PartPrivate.h"
#include "H5Block.h"
#include "H5BlockReadWrite.h"
#include "Underscore.h"

#if defined(F77_SINGLE_UNDERSCORE)
#define F77NAME(a,b) a
#elif defined(F77_CRAY_UNDERSCORE)
#define F77NAME(a,b) b
#elif defined(F77_NO_UNDERSCORE)
#else
#error Error, no way to determine how to construct fortran bindings
#endif
"""

write_scalar_h = """
h5part_int64_t
H5Block#DIM#dWriteScalarField#TYPE_ABV# (
	H5PartFile *f,
	const char *name,
	const h5part_#TYPE_H5P#_t *data
	);
"""

read_scalar_h = """
h5part_int64_t
H5Block#DIM#dReadScalarField#TYPE_ABV# (
	H5PartFile *f,
	const char *name,
	h5part_#TYPE_H5P#_t *data
	);
"""

write_scalar_c = """
/*!
  \\ingroup h5block_data

  Write a 3-dimensional field \\c name from the buffer starting at \\c data
  to the current time-step using the defined field layout. Values are
  #TYPE_FULL#.

  You must use the Fortran indexing scheme to access items in \\c data.

  \\return \\c H5PART_SUCCESS or error code
*/
h5part_int64_t
H5Block#DIM#dWriteScalarField#TYPE_ABV# (
	H5PartFile *f,		  /*!< IN: file handle */
	const char *name,	       /*!< IN: name of dataset to write */
	const h5part_#TYPE_H5P#_t *data      /*!< IN: scalar data to write */
	) {

	SET_FNAME ( "H5Block#DIM#dWriteScalarField#TYPE_ABV#" );
	BLOCK_INIT ( f );
	CHECK_WRITABLE_MODE ( f );
	CHECK_TIMEGROUP ( f );
	CHECK_LAYOUT ( f );

	h5part_int64_t herr = _H5Block_create_field_group ( f, name );
	if ( herr < 0 ) return herr;

	herr = _H5Block_write_data ( f, "0", data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;

	herr = _H5Block_close_field_group ( f );
	if ( herr < 0 ) return herr;

	return H5PART_SUCCESS;
}
"""

read_scalar_c = """
/*!
  \\ingroup h5block_data

  Read a 3-dimensional field \\c name into the buffer starting at \\c data from
  the current time-step using the defined field layout. Values are
  #TYPE_FULL#.

  You must use the Fortran indexing scheme to access items in \\c data.

  \\return \\c H5PART_SUCCESS or error code
*/
h5part_int64_t
H5Block#DIM#dReadScalarField#TYPE_ABV# (
	H5PartFile *f,		  /*!< IN: file handle */
	const char *name,	       /*!< IN: name of dataset to read */
	h5part_#TYPE_H5P#_t *data	  /*!< OUT: ptr to read buffer */
	) {

	SET_FNAME ( "H5Block#DIM#dReadScalarField#TYPE_ABV#" );
	BLOCK_INIT ( f );
	CHECK_TIMEGROUP ( f );
	CHECK_LAYOUT ( f );

	h5part_int64_t herr = _H5Block_open_field_group ( f, name );
	if ( herr < 0 ) return herr;

	herr = _H5Block_read_data ( f, "0", data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;

	herr = _H5Block_close_field_group ( f );
	if ( herr < 0 ) return herr;

	return H5PART_SUCCESS;
}
"""

write_scalar_fi = """
!> \\ingroup h5blockf_data
!! See \\ref H5Block#DIM#dWriteScalarField#TYPE_ABV#
!! \\return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_#DIM#d_write_scalar_field_#TYPE_F90_ABV# ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    #TYPE_F90#, INTENT(IN) :: data(*)    !< the array of data
END FUNCTION
"""

read_scalar_fi = """
!> \\ingroup h5blockf_data
!! See \\ref H5Block#DIM#dReadScalarField#TYPE_ABV#
!! \\return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_#DIM#d_read_scalar_field_#TYPE_F90_ABV# ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    #TYPE_F90#, INTENT(OUT) :: data(*)   !< buffer to read the data into
END FUNCTION
"""

write_scalar_fc = """
#if ! defined(F77_NO_UNDERSCORE)
#define h5bl_#DIM#d_write_scalar_field_#TYPE_F90_ABV# F77NAME ( \\
	h5bl_#DIM#d_write_scalar_field_#TYPE_F90_ABV#_, \\
	H5BL_#DIM#D_WRITE_SCALAR_FIELD_#TYPE_F90_ABVC# )
#endif

h5part_int64_t
h5bl_#DIM#d_write_scalar_field_#TYPE_F90_ABV# (
	h5part_int64_t *f,
	const char *field_name,
	const h5part_#TYPE_H5P#_t *data,
	const int l_field_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *field_name2 =  _H5Part_strdupfor2c ( field_name,  l_field_name );

	h5part_int64_t herr = H5Block#DIM#dWriteScalarField#TYPE_ABV# (
		filehandle, field_name2, data );

	free ( field_name2 );
	return herr;
}
"""

read_scalar_fc = """
#if ! defined(F77_NO_UNDERSCORE)
#define h5bl_3d_read_scalar_field_#TYPE_F90_ABV# F77NAME ( \\
	h5bl_3d_read_scalar_field_#TYPE_F90_ABV#_, \\
	H5BL_#DIM#D_READ_SCALAR_FIELD_#TYPE_F90_ABVC# )
#endif

h5part_int64_t
h5bl_#DIM#d_read_scalar_field_#TYPE_F90_ABV# (
	h5part_int64_t *f,
	const char *field_name,
	h5part_#TYPE_H5P#_t *data,
	const int l_field_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *field_name2 =  _H5Part_strdupfor2c ( field_name,  l_field_name );

	h5part_int64_t herr = H5Block#DIM#dReadScalarField#TYPE_ABV# (
		filehandle, field_name2, data );

	free ( field_name2 );
	return herr;
}
"""

write_vector_h = """
h5part_int64_t
H5Block#DIM#dWrite3dVectorField#TYPE_ABV# (
	H5PartFile *f,
	const char *name,
	const h5part_#TYPE_H5P#_t *xval,
	const h5part_#TYPE_H5P#_t *yval,
	const h5part_#TYPE_H5P#_t *zval
	);
"""

read_vector_h = """
h5part_int64_t
H5Block#DIM#dRead3dVectorField#TYPE_ABV# (
	H5PartFile *f,
	const char *name,
	h5part_#TYPE_H5P#_t *xval,
	h5part_#TYPE_H5P#_t *yval,
	h5part_#TYPE_H5P#_t *zval
	);
"""

write_vector_c = """
/*!
  \\ingroup h5block_data
*/
/*!
  Write a 3-dimensional field \\c name with 3-dimensional vectors as values
  from the buffers starting at \\c x_data, \\c y_data and \\c z_data to the
  current time-step using the defined field layout. Values are 3-dimensional
  vectors with #TYPE_FULL# values.

  You must use the Fortran indexing scheme to access items in \\c data.

  \\return \\c H5PART_SUCCESS or error code
*/
h5part_int64_t
H5Block#DIM#dWrite3dVectorField#TYPE_ABV# (
	H5PartFile *f,		  /*!< IN: file handle */
	const char *name,	       /*!< IN: name of dataset to write */
	const h5part_#TYPE_H5P#_t *x_data, /*!< IN: X axis data */
	const h5part_#TYPE_H5P#_t *y_data, /*!< IN: Y axis data */
	const h5part_#TYPE_H5P#_t *z_data  /*!< IN: Z axis data */
	) {

	SET_FNAME ( "H5Block#DIM#dWrite3dVectorField#TYPE_ABV#" );
	BLOCK_INIT ( f );
	CHECK_WRITABLE_MODE ( f );
	CHECK_TIMEGROUP ( f );
	CHECK_LAYOUT ( f );

	h5part_int64_t herr = _H5Block_create_field_group ( f, name );
	if ( herr < 0 ) return herr;

	herr = _H5Block_write_data ( f, "0", x_data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;
	herr = _H5Block_write_data ( f, "1", y_data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;
	herr = _H5Block_write_data ( f, "2", z_data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;

	herr = _H5Block_close_field_group ( f );
	if ( herr < 0 ) return herr;

	return H5PART_SUCCESS;
}
"""

read_vector_c = """
/*!
  \\ingroup h5block_data
*/
/*!
  Read a 3-dimensional field \\c name with 3-dimensional vectors as values
  from the buffers starting at \\c x_data, \\c y_data and \\c z_data to the
  current time-step using the defined field layout. Values are 3-dimensional
  vectors with #TYPE_FULL# values.

  You must use the Fortran indexing scheme to access items in \\c data.

  \\return \\c H5PART_SUCCESS or error code
*/
h5part_int64_t
H5Block#DIM#dRead3dVectorField#TYPE_ABV# (
	H5PartFile *f,		  /*!< IN: file handle */
	const char *name,	       /*!< IN: name of dataset to write */
	h5part_#TYPE_H5P#_t *x_data, /*!< OUT: X axis data */
	h5part_#TYPE_H5P#_t *y_data, /*!< OUT: Y axis data */
	h5part_#TYPE_H5P#_t *z_data  /*!< OUT: Z axis data */
	) {

	SET_FNAME ( "H5Block#DIM#dRead3dVectorField#TYPE_ABV#" );
	BLOCK_INIT ( f );
	CHECK_TIMEGROUP ( f );
	CHECK_LAYOUT ( f );

	h5part_int64_t herr = _H5Block_open_field_group ( f, name );
	if ( herr < 0 ) return herr;

	herr = _H5Block_read_data ( f, "0", x_data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;
	herr = _H5Block_read_data ( f, "1", y_data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;
	herr = _H5Block_read_data ( f, "2", z_data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;

	herr = _H5Block_close_field_group ( f );
	if ( herr < 0 ) return herr;

	return H5PART_SUCCESS;
}
"""

write_vector_fi = """ 
!> \\ingroup h5blockf_data
!! See \\ref H5Block#DIM#dWrite3dVectorField#TYPE_ABV#
!! \\return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_#DIM#d_write_3dvector_field_#TYPE_F90_ABV# ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    #TYPE_F90#, INTENT(IN) :: x(*)       !< the array of x data to write
    #TYPE_F90#, INTENT(IN) :: y(*)       !< the array of y data to write
    #TYPE_F90#, INTENT(IN) :: z(*)       !< the array of z data to write
END FUNCTION
"""

read_vector_fi = """
!> \\ingroup h5blockf_data
!! See \\ref H5Block#DIM#dRead3dVectorField#TYPE_ABV#
!! \\return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_#DIM#d_read_3dvector_field_#TYPE_F90_ABV# ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    #TYPE_F90#, INTENT(OUT) :: x(*)      !< buffer to read the x data into
    #TYPE_F90#, INTENT(OUT) :: y(*)      !< buffer to read the y data into
    #TYPE_F90#, INTENT(OUT) :: z(*)      !< buffer to read the z data into
END FUNCTION
"""

write_vector_fc = """
#if ! defined(F77_NO_UNDERSCORE)
#define h5bl_#DIM#d_write_3dvector_field_#TYPE_F90_ABV# F77NAME ( \\
	h5bl_#DIM#d_write_3dvector_field_#TYPE_F90_ABV#_, \\
	H5BL_#DIM#D_WRITE_3DVECTOR_FIELD_#TYPE_F90_ABVC# )
#endif

h5part_int64_t
h5bl_#DIM#d_write_3dvector_field_#TYPE_F90_ABV# (
	h5part_int64_t *f,	      /*!< file handle */
	const char *field_name,	 /*!< name of the data set */
	const h5part_#TYPE_H5P#_t *xval,   /*!< array of x component data */
	const h5part_#TYPE_H5P#_t *yval,   /*!< array of y component data */
	const h5part_#TYPE_H5P#_t *zval,   /*!< array of z component data */
	const int l_field_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *field_name2 =  _H5Part_strdupfor2c ( field_name,  l_field_name );

	h5part_int64_t herr = H5Block#DIM#dWrite3dVectorField#TYPE_ABV# (
		filehandle, field_name2, xval, yval, zval );

	free ( field_name2 );
	return herr;
}
"""

read_vector_fc = """
#if ! defined(F77_NO_UNDERSCORE)
#define h5bl_#DIM#d_read_3dvector_field_#TYPE_F90_ABV# F77NAME ( \\
	h5bl_#DIM#d_read_3dvector_field_#TYPE_F90_ABV#_, \\
	H5BL_#DIM#D_READ_3DVECTOR_FIELD_#TYPE_F90_ABVC# )
#endif

h5part_int64_t
h5bl_#DIM#d_read_3dvector_field_#TYPE_F90_ABV# (
	h5part_int64_t *f,	      /*!< file handle */
	const char *field_name,	 /*!< name of the data set */
	h5part_#TYPE_H5P#_t *xval,	 /*!< array of x component data */
	h5part_#TYPE_H5P#_t *yval,	 /*!< array of y component data */
	h5part_#TYPE_H5P#_t *zval,	 /*!< array of z component data */
	const int l_field_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *field_name2 =  _H5Part_strdupfor2c ( field_name,  l_field_name );

	h5part_int64_t herr = H5Block#DIM#dRead3dVectorField#TYPE_ABV# (
		filehandle, field_name2, xval, yval, zval );

	free ( field_name2 );
	return herr;
}
"""

write_attr_h = """
h5part_int64_t
H5BlockWriteFieldAttrib#TYPE_ABV# (
	H5PartFile *f,				/*!< IN: file handle */
	const char *field_name,			/*!< IN: field name */
	const char *attrib_name,		/*!< IN: attribute name */
	const h5part_#TYPE_H5P#_t *attrib_value,		/*!< IN: attribute value */
	const h5part_int64_t attrib_nelem	/*!< IN: number of elements */
	);
"""

write_attr_c = """
/*!
  \\ingroup h5block_attrib

  Write \\c attrib_value with type #TYPE_FULL# as attribute \\c attrib_name
  to field \\c field_name.

  \\return \\c H5PART_SUCCESS or error code
*/
h5part_int64_t
H5BlockWriteFieldAttrib#TYPE_ABV# (
	H5PartFile *f,				/*!< IN: file handle */
	const char *field_name,			/*!< IN: field name */
	const char *attrib_name,		/*!< IN: attribute name */
	const h5part_#TYPE_H5P#_t *attrib_value,		/*!< IN: attribute value */
	const h5part_int64_t attrib_nelem	/*!< IN: number of elements */
	) {

	SET_FNAME ( "H5BlockWriteFieldAttrib#TYPE_ABV#" );
	BLOCK_INIT ( f );
	CHECK_WRITABLE_MODE( f );
	CHECK_TIMEGROUP( f );

	return _write_field_attrib (
		f,
		field_name,
		attrib_name,
                #TYPE_HDF5#,
                attrib_value,
		attrib_nelem );
}
"""

write_attr_fi = """
!> \\ingroup h5blockf_attrib
!! See \\ref H5BlockWriteFieldAttrib#TYPE_ABV#
!! \\return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_writefieldattrib_#TYPE_F90_ABV# ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: field_name  !< the name of the field
    CHARACTER(LEN=*), INTENT(IN) :: attrib_name !< the name of the attribute
    #TYPE_F90#, INTENT(IN) :: attrib_value(*)   !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: attrib_nelem       !< the number of elements in the array
END FUNCTION
"""

write_attr_fc = """
#if ! defined(F77_NO_UNDERSCORE)
#define h5bl_writefieldattrib_#TYPE_F90_ABV# F77NAME ( \\
	h5bl_writefieldattrib_#TYPE_F90_ABV#_, \\
	H5BL_WRITEFIELDATTRIB_#TYPE_F90_ABVC# )
#endif

h5part_int64_t
h5bl_writefieldattrib_#TYPE_F90_ABV# (
	h5part_int64_t *f,
	const char *field_name,
	const char *attrib_name,
	const h5part_#TYPE_H5P#_t *attrib_value,
	const h5part_int64_t *attrib_nelem,
	const int l_field_name,
	const int l_attrib_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *field_name2 = _H5Part_strdupfor2c ( field_name,  l_field_name );
	char *attrib_name2 =_H5Part_strdupfor2c ( attrib_name, l_attrib_name );

	h5part_int64_t herr = H5BlockWriteFieldAttrib#TYPE_ABV# (
		filehandle, field_name2, attrib_name2,
		attrib_value, *attrib_nelem );

	free ( field_name2 );
	free ( attrib_name2 );
	return herr;
}
"""


dims = ["3"]
types = [
  ["floating points (64-bit)", "Float64", "float64", "H5T_NATIVE_DOUBLE", "REAL*8", "r8", "R8"],
  ["floating points (32-bit)", "Float32", "float32", "H5T_NATIVE_FLOAT", "REAL*4", "r4", "R4"],
  ["integers (64-bit)", "Int64", "int64", "H5T_NATIVE_INT64", "INTEGER*8", "i8", "I8"],
  ["integers (32-bit)", "Int32", "int32", "H5T_NATIVE_INT32", "INTEGER*4", "i4", "I4"]
]

def create_call(template, type, dim):
  fcn = template
  fcn = fcn.replace('#DIM#',dim)\
           .replace('#TYPE_FULL#',type[0])\
           .replace('#TYPE_ABV#',type[1])\
           .replace('#TYPE_H5P#',type[2])\
           .replace('#TYPE_HDF5#',type[3])\
           .replace('#TYPE_F90#',type[4])\
           .replace('#TYPE_F90_ABV#',type[5])\
           .replace('#TYPE_F90_ABVC#',type[6]) 
  return fcn

def write_calls():
  cfile = file('H5BlockReadWrite.c','w')
  cfile.write(c_head)
  hfile = file('H5BlockReadWrite.h','w')
  hfile.write(h_head)
  fcfile = file('H5BlockReadWriteF.c','w')
  fcfile.write(fc_head)
  fifile = file('H5BlockReadWrite.f90','w')
  for dim in dims:
    for type in types:
      cfile.write(create_call(write_scalar_c,type,dim));
      cfile.write(create_call(read_scalar_c,type,dim));
      hfile.write(create_call(write_scalar_h,type,dim));
      hfile.write(create_call(read_scalar_h,type,dim));
      fcfile.write(create_call(write_scalar_fc,type,dim));
      fcfile.write(create_call(read_scalar_fc,type,dim));
      fifile.write(create_call(write_scalar_fi,type,dim));
      fifile.write(create_call(read_scalar_fi,type,dim));
      cfile.write(create_call(write_vector_c,type,dim));
      cfile.write(create_call(read_vector_c,type,dim));
      hfile.write(create_call(write_vector_h,type,dim));
      hfile.write(create_call(read_vector_h,type,dim));
      fcfile.write(create_call(write_vector_fc,type,dim));
      fcfile.write(create_call(read_vector_fc,type,dim));
      fifile.write(create_call(write_vector_fi,type,dim));
      fifile.write(create_call(read_vector_fi,type,dim));
  for type in types:
    cfile.write(create_call(write_attr_c,type,""));
    hfile.write(create_call(write_attr_h,type,""));
    fifile.write(create_call(write_attr_fi,type,""));
    fcfile.write(create_call(write_attr_fc,type,""));
  cfile.close()
  hfile.write(h_tail)
  hfile.close()
  fcfile.close()
  fifile.close()

write_calls()

