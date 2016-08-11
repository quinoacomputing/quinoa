#ifndef _H5Part_H_
#define _H5Part_H_

#include <stdlib.h>
#include <stdarg.h>
#include <hdf5.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "H5PartTypes.h"
#include "H5PartAttrib.h"
#include "H5Block.h"
#ifdef PARALLEL_IO
#include "H5MultiBlock.h"
#endif

#define H5PART_VER_STRING	"1.6.6"
#define H5PART_VER_MAJOR	1
#define H5PART_VER_MINOR	6
#define H5PART_VER_RELEASE	6

/* error values */
#define H5PART_SUCCESS		0
#define H5PART_ERR_NOMEM	-12
#define H5PART_ERR_INVAL	-22
#define H5PART_ERR_BADFD	-77

#define H5PART_ERR_INIT         -200
#define H5PART_ERR_NOENTRY	-201
#define H5PART_ERR_NOTYPE	-210
#define H5PART_ERR_BAD_VIEW	-220

#define H5PART_ERR_MPI		-300
#define H5PART_ERR_HDF5		-400

/* file open flags */
#define H5PART_READ		0x01
#define H5PART_WRITE		0x02
#define H5PART_APPEND		0x04
#define H5PART_VFD_MPIPOSIX	0x08
#define H5PART_FS_LUSTRE	0x10
#define H5PART_VFD_MPIIO_IND	0x20
#define H5PART_VFD_CORE		0x40

/* verbosity level flags */
#define H5PART_VERB_NONE	0
#define H5PART_VERB_ERROR	1
#define H5PART_VERB_WARN	2
#define H5PART_VERB_INFO	3
#define H5PART_VERB_DEBUG	4
#define H5PART_VERB_DETAIL	5

/* data types */
#define H5PART_INT64		((h5part_int64_t)H5T_NATIVE_INT64)
#define H5PART_INT32		((h5part_int64_t)H5T_NATIVE_INT32) 
#define H5PART_FLOAT64		((h5part_int64_t)H5T_NATIVE_DOUBLE)
#define H5PART_FLOAT32		((h5part_int64_t)H5T_NATIVE_FLOAT)
#define H5PART_CHAR		((h5part_int64_t)H5T_NATIVE_CHAR)
#define H5PART_STRING		((h5part_int64_t)H5T_C_S1)

/*========== File Opening/Closing ===============*/
H5PartFile*
H5PartOpenFile(
	const char *filename,
	const char flags
	);

H5PartFile*
H5PartOpenFileAlign(
	const char *filename,
	const char flags,
	h5part_int64_t align
	);

#define H5PartOpenFileSerial(x,y) H5PartOpenFile(x,y)
#define H5PartOpenFileSerialAlign(x,y,z) H5PartOpenFileAlign(x,y,z)

#ifdef PARALLEL_IO
H5PartFile*
H5PartOpenFileParallel (
	const char *filename,
	const char flags,
	MPI_Comm communicator
	);

H5PartFile*
H5PartOpenFileParallelAlign (
	const char *filename,
	const char flags,
	MPI_Comm communicator,
	h5part_int64_t align
	);
#endif


h5part_int64_t
H5PartCloseFile (
	H5PartFile *f
	);

h5part_int64_t
H5PartFileIsValid (
	H5PartFile *f
	);

/*============== File Writing Functions ==================== */
h5part_int64_t
H5PartDefineStepName (
	H5PartFile *f,
	const char *name,
	const h5part_int64_t width
	);

h5part_int64_t
H5PartSetNumParticles ( 
	H5PartFile *f, 
	const h5part_int64_t nparticles
	);

h5part_int64_t
H5PartSetNumParticlesStrided (
	H5PartFile *f,				/*!< [in] Handle to open file */
	const h5part_int64_t nparticles,	/*!< [in] Number of particles */
	const h5part_int64_t stride		/*!< [in] Stride (e.g. number of fields in the particle array) */
	);

h5part_int64_t
H5PartSetChunkSize (
	H5PartFile *f,
	const h5part_int64_t size
	);

h5part_int64_t
H5PartWriteDataFloat64 (
	H5PartFile *f,
	const char *name,
	const h5part_float64_t *array
	);

h5part_int64_t
H5PartWriteDataFloat32 (
	H5PartFile *f,
	const char *name,
	const h5part_float32_t *array
	);

h5part_int64_t
H5PartWriteDataInt64 (
	H5PartFile *f,
	const char *name,
	const h5part_int64_t *array
	);

h5part_int64_t
H5PartWriteDataInt32 (
	H5PartFile *f,
	const char *name,
	const h5part_int32_t *array
	);

/*================== File Reading Routines =================*/
h5part_int64_t
H5PartSetStep (
	H5PartFile *f,
	const h5part_int64_t step
	);

h5part_int64_t
H5PartHasStep (
	H5PartFile *f,
	const h5part_int64_t step
	);

h5part_int64_t
H5PartGetNumSteps (
	H5PartFile *f
	);

h5part_int64_t
H5PartGetNumDatasets (
	H5PartFile *f
	);

h5part_int64_t
H5PartGetDatasetName (
	H5PartFile *f,
	const h5part_int64_t idx,
	char *name,
	const h5part_int64_t maxlen
	);

h5part_int64_t
H5PartGetDatasetInfo (
	H5PartFile *f,
	const h5part_int64_t idx,
	char *name,
	const h5part_int64_t maxlen,
	h5part_int64_t *type,
	h5part_int64_t *nelem);


h5part_int64_t
H5PartGetNumParticles (
	H5PartFile *f
	);

h5part_int64_t
H5PartSetView (
	H5PartFile *f,
	const h5part_int64_t start,
	const h5part_int64_t end
	);

h5part_int64_t
H5PartSetViewIndices (
	H5PartFile *f,			/*!< [in]  Handle to open file */
	const h5part_int64_t *indices,	/*!< [in]  List of indices */
	h5part_int64_t nelems		/*!< [in]  Size of list */
	);

h5part_int64_t
H5PartGetView (
	H5PartFile *f,
	h5part_int64_t *start,
	h5part_int64_t *end
	);

h5part_int64_t
H5PartHasView (
	H5PartFile *f
	);

h5part_int64_t
H5PartResetView (
	H5PartFile *f
	);

h5part_int64_t
H5PartSetCanonicalView (
	H5PartFile *f
	);

h5part_int64_t
H5PartReadDataFloat64(
	H5PartFile *f,
	const char *name,
	h5part_float64_t *array
	);

h5part_int64_t
H5PartReadDataFloat32(
	H5PartFile *f,
	const char *name,
	h5part_float32_t *array
	);

h5part_int64_t
H5PartReadDataInt64 (
	H5PartFile *f,
	const char *name,
	h5part_int64_t *array
	);

h5part_int64_t
H5PartReadDataInt32 (
	H5PartFile *f,
	const char *name,
	h5part_int32_t *array
	);

h5part_int64_t
H5PartReadParticleStep (
	H5PartFile *f,
	const h5part_int64_t step,
	h5part_float64_t *x, /* particle positions */
	h5part_float64_t *y,
	h5part_float64_t *z,
	h5part_float64_t *px, /* particle momenta */
	h5part_float64_t *py,
	h5part_float64_t *pz,
	h5part_int64_t *id /* and phase */
	);

/**********==============Attributes Interface============***************/
/* currently there is file attributes:  Attributes bound to the file
   and step attributes which are bound to the current timestep.  You 
   must set the timestep explicitly before writing the attributes (just
   as you must do when you write a new dataset.  Currently there are no
   attributes that are bound to a particular data array, but this could
   easily be done if required.
*/
h5part_int64_t
H5PartWriteStepAttrib (
	H5PartFile *f,
	const char *attrib_name,
	const h5part_int64_t attrib_type,
	const void *attrib_value,
	const h5part_int64_t attrib_nelem
	);

h5part_int64_t
H5PartWriteFileAttrib (
	H5PartFile *f,
	const char *attrib_name,
	const h5part_int64_t attrib_type,
	const void *attrib_value,
	const h5part_int64_t attrib_nelem
	);

h5part_int64_t
H5PartWriteFileAttribString (
	H5PartFile *f,
	const char *name,
	const char *value
	);

h5part_int64_t
H5PartWriteStepAttribString ( 
	H5PartFile *f,
	const char *name,
	const char *value
	);

h5part_int64_t
H5PartGetNumStepAttribs ( /* for current filestep */
	H5PartFile *f
	);

h5part_int64_t
H5PartGetNumFileAttribs (
	H5PartFile *f
	);

h5part_int64_t
H5PartGetStepAttribInfo (
	H5PartFile *f,
	const h5part_int64_t attrib_idx,
	char *attrib_name,
	const h5part_int64_t len_of_attrib_name,
	h5part_int64_t *attrib_type,
	h5part_int64_t *attrib_nelem
	);

h5part_int64_t
H5PartGetFileAttribInfo (
	H5PartFile *f,
	const h5part_int64_t idx,
	char *name,
	const h5part_int64_t maxnamelen,
	h5part_int64_t *type,
	h5part_int64_t *nelem
	);

h5part_int64_t
H5PartReadStepAttrib (
	H5PartFile *f,
	const char *name,
	void *value
	);

h5part_int64_t
H5PartReadFileAttrib (
	H5PartFile *f,
	const char *name,
	void *value
	);

/*============ Error Reporting and Configuration =============*/

h5part_int64_t
H5PartSetVerbosityLevel (
	const h5part_int64_t level
	);

h5part_int64_t
H5PartSetThrottle (
	H5PartFile *f,
	int factor
	);

h5part_int64_t
H5PartSetErrorHandler (
	const h5part_error_handler handler
	);

h5part_int64_t
H5PartGetErrno (
	void
	);

h5part_error_handler
H5PartGetErrorHandler (
	void
	);

h5part_int64_t
H5PartReportErrorHandler (
	const char *funcname,
	const h5part_int64_t eno,
	const char *fmt,
	...
	);

h5part_int64_t
H5PartAbortErrorHandler (
	const char *funcname,
	const h5part_int64_t eno,
	const char *fmt,
	...
	);

#ifdef __cplusplus
}
#endif

#endif
