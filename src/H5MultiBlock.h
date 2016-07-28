#ifndef _H5MULTIBLOCK_H_
#define _H5MULTIBLOCK_H_

#ifdef __cplusplus
extern "C" {
#endif

/*!
  Include read/write call variants for different data types and
  field dimensions.
*/
#include "H5MultiBlockReadWrite.h"

#define H5MULTIBLOCK_ATTR_NAME "__BlockDims__"

h5part_int64_t
H5MultiBlock3dDefineRadius (
	H5PartFile *f,
	const h5part_int64_t r
	);

h5part_int64_t
H5MultiBlock3dDefineRadii (
	H5PartFile *f,
	const h5part_int64_t ri,
	const h5part_int64_t rj,
	const h5part_int64_t rk
	);

h5part_int64_t
H5MultiBlock3dDefineDims (
	H5PartFile *f,
	const h5part_int64_t *field_dims,
	const h5part_int64_t *block_dims
	);

h5part_int64_t
H5MultiBlock3dGetFieldDims(
	H5PartFile *f,
	h5part_int64_t *dims
	);

h5part_int64_t
H5MultiBlock3dGetBlockDims(
	H5PartFile *f,
	const char *field_name,
	h5part_int64_t *dims
	);

h5part_int64_t
H5MultiBlock3dGetOffsetsOfProc (
	H5PartFile *f,
	const h5part_int64_t proc,
	h5part_int64_t *offsets
	);

h5part_int64_t
H5MultiBlock3dCalculateDecomp (
	const int nprocs,
	h5part_int64_t *decomp
	);

h5part_int64_t
H5MultiBlockFree (
	void *block
	);

h5part_int64_t
H5MultiBlockShiftProcs (
	H5PartFile *f,
	const int shift
	);

#ifdef __cplusplus
}
#endif

#endif
