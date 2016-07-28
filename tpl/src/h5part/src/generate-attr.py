#!/usr/bin/python

c_head = """
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>
#include "H5Part.h"
#include "H5PartErrors.h"
#include "H5PartPrivate.h"

"""

h_head = """
#ifndef _H5PART_ATTRIB_H_
#define _H5PART_ATTRIB_H_

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

write_attr_h = """
h5part_int64_t
H5PartWrite#LEVEL#Attrib#TYPE_ABV# (
	H5PartFile *f,
	const char *name,
	const h5part_#TYPE_H5P#_t data
	);
"""

write_attr_c = """
/*!
  \\ingroup h5part_attrib

  Writes a \\c value of type #TYPE_FULL#
  to the root ("/") of the file
  as attribute \\c name.

  \\return \\c H5PART_SUCCESS or error code
*/
h5part_int64_t
H5PartWrite#LEVEL#Attrib#TYPE_ABV# (
	H5PartFile *f,				/*!< IN: file handle */
	const char *name,		        /*!< IN: attribute name */
	const h5part_#TYPE_H5P#_t value	        /*!< IN: attribute value */
	) {

	SET_FNAME ( "H5PartWrite#LEVEL#Attrib#TYPE_ABV#" );

	CHECK_FILEHANDLE ( f );
	CHECK_WRITABLE_MODE( f );

	h5part_int64_t herr = _H5Part_write_#LEVELLC#_attrib (
		f,
		name,
                #TYPE_HDF5#,
                (void*)&value,
		1 );
	if ( herr < 0 ) return herr;

	return H5PART_SUCCESS;
}
"""

write_attr_fi = """
!< \\ingroup h5partf_attrib
!! See \\ref H5PartWrite#LEVELLC#Attrib#TYPE_ABV#
!! \\return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_write#LEVELLC#attrib_#TYPE_F90_ABV# ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    #TYPE_F90#, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION
"""

write_attr_fc = """
#if ! defined(F77_NO_UNDERSCORE)
#define h5pt_write#LEVELLC#attrib_#TYPE_F90_ABV# F77NAME ( \\
	h5pt_write#LEVELLC#attrib_#TYPE_F90_ABV#_, \\
	H5PT_WRITE#LEVELUC#ATTRIB_#TYPE_F90_ABVC# )
#endif

h5part_int64_t
h5pt_write#LEVELLC#attrib_#TYPE_F90_ABV# (
	h5part_int64_t *f,
	const char *name,
	const h5part_#TYPE_H5P#_t *data,
	const h5part_int64_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWrite#LEVEL#Attrib (
		filehandle, name2, #TYPE_ATTRIB#, data, *nelem);

	free ( name2 );
	return herr;
}
"""

read_attr_fi = """
!< \\ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \\return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_read#LEVELLC#attrib_#TYPE_F90_ABV# ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    #TYPE_F90#, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION
"""

read_attr_fc = """
#if ! defined(F77_NO_UNDERSCORE)
#define h5pt_read#LEVELLC#attrib_#TYPE_F90_ABV# F77NAME ( \\
	h5pt_read#LEVELLC#attrib_#TYPE_F90_ABV#_, \\
	H5PT_READ#LEVELUC#ATTRIB_#TYPE_F90_ABVC# )
#endif

h5part_int64_t
h5pt_read#LEVELLC#attrib_#TYPE_F90_ABV# (
	h5part_int64_t *f,
	const char *name,
	const h5part_#TYPE_H5P#_t *data,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartRead#LEVEL#Attrib (
		filehandle, name2, (void*)data);

	free ( name2 );
	return herr;
}
"""


levels = ["File", "Step"]
types = [
  ["floating points (64-bit)", "Float64", "float64", "H5T_NATIVE_DOUBLE", "REAL*8", "r8", "R8", "H5PART_FLOAT64"],
  ["floating points (32-bit)", "Float32", "float32", "H5T_NATIVE_FLOAT", "REAL*4", "r4", "R4", "H5PART_FLOAT32"],
  ["integers (64-bit)", "Int64", "int64", "H5T_NATIVE_INT64", "INTEGER*8", "i8", "I8", "H5PART_INT64"],
  ["integers (32-bit)", "Int32", "int32", "H5T_NATIVE_INT32", "INTEGER*4", "i4", "I4", "H5PART_INT32"]
]

def create_call(template, type, level):
  fcn = template
  fcn = fcn.replace('#LEVEL#',level)\
           .replace('#LEVELLC#',level.lower())\
           .replace('#LEVELUC#',level.upper())\
           .replace('#TYPE_FULL#',type[0])\
           .replace('#TYPE_ABV#',type[1])\
           .replace('#TYPE_H5P#',type[2])\
           .replace('#TYPE_HDF5#',type[3])\
           .replace('#TYPE_F90#',type[4])\
           .replace('#TYPE_F90_ABV#',type[5])\
           .replace('#TYPE_F90_ABVC#',type[6])\
           .replace('#TYPE_ATTRIB#',type[7]) 
  return fcn

def write_calls():
  cfile = file('H5PartAttrib.c','w')
  cfile.write(c_head)
  hfile = file('H5PartAttrib.h','w')
  hfile.write(h_head)
  fcfile = file('H5PartAttribF.c','w')
  fcfile.write(fc_head)
  fifile = file('H5PartAttrib.f90','w')
  for level in levels:
    for type in types:
      cfile.write(create_call(write_attr_c,type,level));
      hfile.write(create_call(write_attr_h,type,level));
      fifile.write(create_call(write_attr_fi,type,level));
      fcfile.write(create_call(write_attr_fc,type,level));
      fifile.write(create_call(read_attr_fi,type,level));
      fcfile.write(create_call(read_attr_fc,type,level));
  cfile.close()
  hfile.write(h_tail)
  hfile.close()
  fcfile.close()
  fifile.close()

write_calls()

