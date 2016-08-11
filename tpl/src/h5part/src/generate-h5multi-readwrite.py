#!/usr/bin/python

c_head = """
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>
#include "H5Part.h"
#include "H5PartErrors.h"
#include "H5PartPrivate.h"

#include "H5MultiBlockErrors.h"
#include "H5MultiBlockPrivate.h"

#ifdef PARALLEL_IO

"""

c_tail = """
#endif // PARALLEL_IO
"""

h_head = """
#ifndef _H5MULTIBLOCK_READWRITE_H_
#define _H5MULTIBLOCK_READWRITE_H_

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

write_h = """
h5part_int64_t
H5MultiBlock#DIM#dWriteField#TYPE_ABV# (
	H5PartFile *f,
	const char *name,
	const h5part_#TYPE_H5P#_t *data
	);
"""

read_h = """
h5part_int64_t
H5MultiBlock#DIM#dReadField#TYPE_ABV# (
	H5PartFile *f,
	const char *name,
	h5part_#TYPE_H5P#_t **data
	);
"""

write_c = """
/*!
  \\ingroup h5multiblock_data

  Write a multiblock field \\c name from the buffer starting at \\c data
  to the current time-step using the defined block decomposition and dimensions.
  Values are #TYPE_FULL#.

  You must use the Fortran indexing scheme to access items in \\c data.

  \\return \\c H5PART_SUCCESS or error code
*/
h5part_int64_t
H5MultiBlock#DIM#dWriteField#TYPE_ABV# (
	H5PartFile *f,		/*!< IN: file handle */
	const char *name,       /*!< IN: name of dataset to write */
	const h5part_#TYPE_H5P#_t *data      /*!< IN: data to write */
	) {

	SET_FNAME( "H5MultiBlock#DIM#dWriteField#TYPE_ABV#" );

	h5part_int64_t herr;

	herr = _H5MultiBlock_write_data ( f, name, data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;

	return H5PART_SUCCESS;
}
"""

read_c = """
/*!
  \\ingroup h5multiblock_data

  Allocate a buffer to hold a block from a multiblock field and place the
  pointer in \\c data, then read the block into the buffer. Uses the block
  decomposition specified in the file and the defined halo radius.
  Values are #TYPE_FULL#.

  You must use the Fortran indexing scheme to access items in \\c data.

  \\return \\c H5PART_SUCCESS or error code
*/
h5part_int64_t
H5MultiBlock#DIM#dReadField#TYPE_ABV# (
	H5PartFile *f,		/*!< IN: file handle */
	const char *name,	/*!< IN: name of dataset to read */
	h5part_#TYPE_H5P#_t **data	/*!< OUT: ptr to read buffer */
	) {

	SET_FNAME( "H5MultiBlock#DIM#dReadField#TYPE_ABV#" );

	h5part_int64_t herr;

	herr = _H5MultiBlock_read_data ( f, name, (char**) data, #TYPE_HDF5# );
	if ( herr < 0 ) return herr;

	return H5PART_SUCCESS;
}
"""

write_fi = """
INTEGER*8 FUNCTION h5multi_#DIM#d_write_field_#TYPE_F90_ABV# ( filehandle, name, data )
  INTEGER*8, INTENT(IN) :: filehandle
  CHARACTER(LEN=*), INTENT(IN) :: name
  #TYPE_F90#, INTENT(IN) :: data(*)
END FUNCTION
"""

read_fi = """
INTEGER*8 FUNCTION h5multi_#DIM#d_read_field_#TYPE_F90_ABV# ( filehandle, name, data )
  INTEGER*8, INTENT(IN) :: filehandle
  CHARACTER(LEN=*), INTENT(IN) :: name
  #TYPE_F90#, INTENT(OUT) :: data(*)
END FUNCTION
"""

write_fc = """
#if ! defined(F77_NO_UNDERSCORE)
#define h5multi_#DIM#d_write_field_#TYPE_F90_ABV# F77NAME ( \\
	h5multi_#DIM#d_write_field_#TYPE_F90_ABV#_, \\
	H5MULTI_#DIM#D_WRITE_FIELD_#TYPE_F90_ABVC# )
#endif

h5part_int64_t
h5multi_#DIM#d_write_field_#TYPE_F90_ABV# (
	h5part_int64_t *f,
	const char *field_name,
	const h5part_#TYPE_H5P#_t *data,
	const int l_field_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *field_name2 = _H5Part_strdupfor2c ( field_name,  l_field_name );

	h5part_int64_t herr = H5MultiBlock#DIM#dWriteField#TYPE_ABV# (
		filehandle, field_name2, data );

	free ( field_name2 );
	return herr;
}
"""

read_fc = """
#if ! defined(F77_NO_UNDERSCORE)
#define h5multi_#DIM#d_read_field_#TYPE_F90_ABV# F77NAME ( \\
	h5multi_#DIM#d_read_field_#TYPE_F90_ABV#_, \\
	H5MULTI_#DIM#D_READ_FIELD_#TYPE_F90_ABVC# )
#endif

h5part_int64_t
h5multi_#DIM#d_read_field_#TYPE_F90_ABV# (
	h5part_int64_t *f,
	const char *field_name,
	h5part_#TYPE_H5P#_t **data,
	const int l_field_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *field_name2 = _H5Part_strdupfor2c ( field_name,  l_field_name );

	h5part_int64_t herr = H5MultiBlock#DIM#dReadField#TYPE_ABV# (
		filehandle, field_name2, (char**) data );

	free ( field_name2 );
	return herr;
}
"""

dims = [ "3" ]
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
  cfile = file('H5MultiBlockReadWrite.c','w')
  cfile.write(c_head)
  hfile = file('H5MultiBlockReadWrite.h','w')
  hfile.write(h_head)
#  fcfile = file('H5MultiBlockReadWriteF.c','w')
#  fcfile.write(fc_head)
#  fifile = file('H5MultiBlockReadWriteF90.inc','w')
#  fifile.write(fi_head)
  for dim in dims:
    for type in types:
	cfile.write(create_call(write_c,type,dim));
	cfile.write(create_call(read_c,type,dim));
	hfile.write(create_call(write_h,type,dim));
	hfile.write(create_call(read_h,type,dim));
#	fcfile.write(create_call(write_fc,type,dim));
#	fcfile.write(create_call(read_fc,type,dim));
#	fifile.write(create_call(write_fi,type,dim));
#	fifile.write(create_call(read_fi,type,dim));
  cfile.write(c_tail)
  cfile.close()
  hfile.write(h_tail)
  hfile.close()
#  fcfile.close()
#  fifile.write(fi_tail)
#  fifile.close()

write_calls()

