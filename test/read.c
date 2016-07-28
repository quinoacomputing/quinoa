#include <string.h>
#include <H5Part.h>

#include "testframe.h"
#include "params.h"

static void
test_read_file_attribs(H5PartFile *file, int position)
{
	h5part_int64_t status;
	char name[ATTR_NAME_SIZE];

	char str[ATTR_NAME_SIZE];
	h5part_int32_t i32;
	h5part_int64_t i64;
	h5part_float32_t f32;
	h5part_float64_t f64;

	TEST("Reading file attributes");

	i64 = H5PartGetNumFileAttribs(file);
	VALUE(i64 % 5, 0, "file attribute count");

	get_attr_name(name, "str", position);
	status = H5PartReadFileAttrib(file, name, str);
	RETURN(status, H5PART_SUCCESS, "H5PartReadFileAttrib");
	SVALUE(str, ATTR_STR_VAL, "string attribute");

	get_attr_name(name, "i32", position);
	status = H5PartReadFileAttrib(file, name, &i32);
	RETURN(status, H5PART_SUCCESS, "H5PartReadFileAttrib");
	IVALUE(i32, ATTR_INT32_VAL, "int32 attribute");

	get_attr_name(name, "i64", position);
	status = H5PartReadFileAttrib(file, name, &i64);
	RETURN(status, H5PART_SUCCESS, "H5PartReadFileAttrib");
	IVALUE(i64, ATTR_INT64_VAL, "int64 attribute");

	get_attr_name(name, "f32", position);
	status = H5PartReadFileAttrib(file, name, &f32);
	RETURN(status, H5PART_SUCCESS, "H5PartReadFileAttrib");
	FVALUE(f32, ATTR_FLOAT_VAL, "float32 attribute");

	get_attr_name(name, "f64", position);
	status = H5PartReadFileAttrib(file, name, &f64);
	RETURN(status, H5PART_SUCCESS, "H5PartReadFileAttrib");
	FVALUE(f64, ATTR_FLOAT_VAL, "float64 attribute");
}

static void
test_read_step_attribs(H5PartFile *file, int position)
{
	h5part_int64_t status;
	char name[ATTR_NAME_SIZE];

	char str[ATTR_NAME_SIZE];
	h5part_int32_t i32;
	h5part_int64_t i64;
	h5part_float32_t f32;
	h5part_float64_t f64;

	TEST("Reading step attributes");

	i64 = H5PartGetNumStepAttribs(file);
	IVALUE(i64, 5, "step attribute count");

	get_attr_name(name, "str", position);
	status = H5PartReadStepAttrib(file, name, str);
	RETURN(status, H5PART_SUCCESS, "H5PartReadStepAttrib");
	SVALUE(str, ATTR_STR_VAL, "string attribute");

	get_attr_name(name, "i32", position);
	status = H5PartReadStepAttrib(file, name, &i32);
	RETURN(status, H5PART_SUCCESS, "H5PartReadStepAttrib");
	IVALUE(i32, ATTR_INT32_VAL, "int32 attribute");

	get_attr_name(name, "i64", position);
	status = H5PartReadStepAttrib(file, name, &i64);
	RETURN(status, H5PART_SUCCESS, "H5PartReadStepAttrib");
	IVALUE(i64, ATTR_INT64_VAL, "int64 attribute");

	get_attr_name(name, "f32", position);
	status = H5PartReadStepAttrib(file, name, &f32);
	RETURN(status, H5PART_SUCCESS, "H5PartReadStepAttrib");
	FVALUE(f32, ATTR_FLOAT_VAL, "float32 attribute");

	get_attr_name(name, "f64", position);
	status = H5PartReadStepAttrib(file, name, &f64);
	RETURN(status, H5PART_SUCCESS, "H5PartReadStepAttrib");
	FVALUE(f64, ATTR_FLOAT_VAL, "float64 attribute");
}

static void
test_read_data64(H5PartFile *file, int nparticles, int step)
{
	int i,t;
	int rank, nprocs;
	h5part_int64_t status, val, start, end, type, size;
	char name1[4];
	char name2[8];
	h5part_int64_t indices[8];

	double *x,*y,*z;
	double *px,*py,*pz;
	h5part_int64_t *id;

	x=(double*)malloc(nparticles*sizeof(double));
	y=(double*)malloc(nparticles*sizeof(double));
	z=(double*)malloc(nparticles*sizeof(double));
	px=(double*)malloc(nparticles*sizeof(double));
	py=(double*)malloc(nparticles*sizeof(double));
	pz=(double*)malloc(nparticles*sizeof(double));
	id=(h5part_int64_t*)malloc(nparticles*sizeof(h5part_int64_t));

	TEST("Verifying dataset info");

#if PARALLEL_IO
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#else
	nprocs = 1;
	rank = 2;
#endif

	val = H5PartGetNumParticles(file);
	IVALUE(val, nprocs*nparticles, "particle count");

	val = H5PartGetNumDatasets(file);
	IVALUE(val, 7, "dataset count");

	for (i=0; i<7; i++) {
	        status = H5PartGetDatasetName(file, i, name1, 2);
		RETURN(status, H5PART_SUCCESS, "H5PartGetDatasetName");

		status = H5PartGetDatasetInfo(
			file, i, name2, 4, &type, &size);
		RETURN(status, H5PART_SUCCESS, "H5PartGetDatasetInfo");
		CVALUE(name1[0], name2[0], "dataset name");

	        status = H5PartGetDatasetName(file, i, name1, 4);
		RETURN(status, H5PART_SUCCESS, "H5PartGetDatasetName");
		CVALUE(name1[1], name2[1], "dataset name");

		IVALUE(size, nprocs*nparticles, "dataset size");
		if (name1[0] == 'i') IVALUE(type, H5PART_INT64, "dataset type");
		else IVALUE(type, H5PART_FLOAT64, "dataset type");
	}

#if PARALLEL_IO
	TEST("Setting throttle");
	status = H5PartSetThrottle(file, 3);
	RETURN(status, H5PART_SUCCESS, "H5PartSetThrottle");
#endif

	TEST("Reading 64-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		val = H5PartHasStep(file, t);
		IVALUE(val, 1, "has step");

		status = H5PartSetStep(file, t);
		RETURN(status, H5PART_SUCCESS, "H5PartSetStep");

                test_read_step_attribs(file, t);

		status = H5PartSetNumParticles(file, nparticles);
		RETURN(status, H5PART_SUCCESS, "H5PartSetNumParticles");

		status = H5PartResetView(file);
		RETURN(status, H5PART_SUCCESS, "H5PartResetView");

		start = rank;
		end = -1;

		status = H5PartSetView(file, start, end);
		RETURN(status, H5PART_SUCCESS, "H5PartSetView");

		val = H5PartGetView(file, &start, &end);
		IVALUE(val, nprocs*nparticles-start, "particle count");
		IVALUE(start, rank, "view start");
		IVALUE(end, nprocs*nparticles-1, "view end");

		status = H5PartSetView(file, -1, -1);
		RETURN(status, H5PART_SUCCESS, "H5PartSetView");

		status = H5PartSetView(file, 0, nparticles-1);
		RETURN(status, H5PART_SUCCESS, "H5PartSetView");

		val = H5PartGetNumParticles(file);
		IVALUE(val, nparticles, "particle count");

		status = H5PartReadDataFloat64(file, "x", x);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat64");
		IVALUE(x[rank], (double)(rank+nparticles*t), "x data");

		status = H5PartResetView(file);
		RETURN(status, H5PART_SUCCESS, "H5PartResetView");

		val = H5PartGetNumParticles(file);
		IVALUE(val, nprocs*nparticles, "particle count");

		indices[0] = rank*2 + 0;
		indices[1] = rank*2 + 3;
		indices[2] = rank*2 + 9;
		indices[3] = rank*2 + 7;

		status = H5PartSetViewIndices(file, indices, 4);
		RETURN(status, H5PART_SUCCESS, "H5PartSetViewIndices");

		val = H5PartGetNumParticles(file);
		IVALUE(val, 4, "particle count");

		double x2[4];
		status = H5PartReadDataFloat64(file, "x", x2);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat64");
		FVALUE(x2[0], (double)(2*rank+0+nparticles*t), "x data");
		FVALUE(x2[1], (double)(2*rank+3+nparticles*t), "x data");
		FVALUE(x2[2], (double)(2*rank+9+nparticles*t), "x data");
		FVALUE(x2[3], (double)(2*rank+7+nparticles*t), "x data");

		status = H5PartSetViewIndices(file, indices, -1);
		RETURN(status, H5PART_SUCCESS, "H5PartSetViewIndices");

		val = H5PartGetNumParticles(file);
		IVALUE(val, nprocs*nparticles, "particle count");

		status = H5PartSetCanonicalView(file);
		RETURN(status, H5PART_SUCCESS, "H5PartSetCanonicalView");

		val = H5PartGetNumParticles(file);
		IVALUE(val, nparticles, "particle count");

                status = H5PartReadParticleStep (
                        file, t, x, y, z, px, py, pz, id);
		RETURN(status, H5PART_SUCCESS, "H5PartReadParticleStep");

		status = H5PartSetViewIndices(file, NULL, 4);
		RETURN(status, H5PART_SUCCESS, "H5PartSetViewIndices");

		for (i=0; i<nparticles; i++)
		{
			FVALUE(x[i] , 0.0 + (double)(i+nparticles*t), " x data");
			FVALUE(y[i] , 0.1 + (double)(i+nparticles*t), " y data");
			FVALUE(z[i] , 0.2 + (double)(i+nparticles*t), " z data");
			FVALUE(px[i], 0.3 + (double)(i+nparticles*t), " px data");
			FVALUE(py[i], 0.4 + (double)(i+nparticles*t), " py data");
			FVALUE(pz[i], 0.5 + (double)(i+nparticles*t), " pz data");
			IVALUE(id[i],               (i+nparticles*t), " id data");
		}
	}
}

static void
test_read_strided_data64(H5PartFile *file, int nparticles, int step)
{
	int i,t;
	h5part_int64_t status;

	double *data;

	data=(double*)malloc(6*nparticles*sizeof(double));

	TEST("Reading 64-bit strided data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		status = H5PartSetStep(file, t);
		RETURN(status, H5PART_SUCCESS, "H5PartSetStep");

		status = H5PartSetNumParticlesStrided(file, nparticles, 6);
		RETURN(status, H5PART_SUCCESS, "H5PartSetNumParticlesStrided");

		status = H5PartReadDataFloat64(file, "x", data+0);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "y", data+1);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "z", data+2);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat64");
		
		status = H5PartReadDataFloat64(file, "px", data+3);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat64");

                test_read_step_attribs(file, t);

		status = H5PartReadDataFloat64(file, "py", data+4);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat64");

		status = H5PartReadDataFloat64(file, "pz", data+5);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat64");

		for (i=0; i<nparticles; i++)
		{
			FVALUE(data[6*i]  , 0.0 + (double)(i+nparticles*t), "x data");
			FVALUE(data[6*i+1], 0.1 + (double)(i+nparticles*t), "y data");
			FVALUE(data[6*i+2], 0.2 + (double)(i+nparticles*t), "z data");
			FVALUE(data[6*i+3], 0.3 + (double)(i+nparticles*t), "px data");
			FVALUE(data[6*i+4], 0.4 + (double)(i+nparticles*t), "py data");
			FVALUE(data[6*i+5], 0.5 + (double)(i+nparticles*t), "pz data");
		}
	}
}

static void
test_read_data32(H5PartFile *file, int nparticles, int step)
{
	int i,t;
	h5part_int64_t status, val;

	float *x,*y,*z;
	float *px,*py,*pz;
	int *id;

	x=(float*)malloc(nparticles*sizeof(float));
	y=(float*)malloc(nparticles*sizeof(float));
	z=(float*)malloc(nparticles*sizeof(float));
	px=(float*)malloc(nparticles*sizeof(float));
	py=(float*)malloc(nparticles*sizeof(float));
	pz=(float*)malloc(nparticles*sizeof(float));
	id=(int*)malloc(nparticles*sizeof(int));

	TEST("Reading 32-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		val = H5PartHasStep(file, t);
		IVALUE(val, 1, "has step");

		status = H5PartSetStep(file, t);
		RETURN(status, H5PART_SUCCESS, "H5PartSetStep");

		status = H5PartSetCanonicalView(file);
		RETURN(status, H5PART_SUCCESS, "H5PartSetCanonicalView");

                test_read_step_attribs(file, t);

		status = H5PartReadDataFloat32(file, "x", x);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "y", y);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "z", z);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");
		
		status = H5PartReadDataFloat32(file, "px", px);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "py", py);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "pz", pz);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");
		
		status = H5PartReadDataInt32(file, LONGNAME, id);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataInt32");

		for (i=0; i<nparticles; i++)
		{
			FVALUE(x[i] , 0.0F + (float)(i+nparticles*t), " x data");
			FVALUE(y[i] , 0.1F + (float)(i+nparticles*t), " y data");
			FVALUE(z[i] , 0.2F + (float)(i+nparticles*t), " z data");
			FVALUE(px[i], 0.3F + (float)(i+nparticles*t), " px data");
			FVALUE(py[i], 0.4F + (float)(i+nparticles*t), " py data");
			FVALUE(pz[i], 0.5F + (float)(i+nparticles*t), " pz data");
			IVALUE(id[i],               (i+nparticles*t), " id data");
		}
	}
}

static void
test_read_strided_data32(H5PartFile *file, int nparticles, int step)
{
	int i,t;
	h5part_int64_t status;

	float *data;

	data=(float*)malloc(6*nparticles*sizeof(float));

	TEST("Reading 32-bit strided data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		status = H5PartSetStep(file, t);
		RETURN(status, H5PART_SUCCESS, "H5PartSetStep");

		status = H5PartSetNumParticlesStrided(file, nparticles, 6);
		RETURN(status, H5PART_SUCCESS, "H5PartSetNumParticlesStrided");

		status = H5PartReadDataFloat32(file, "x", data+0);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "y", data+1);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "z", data+2);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");
		
		status = H5PartReadDataFloat32(file, "px", data+3);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "py", data+4);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");

		status = H5PartReadDataFloat32(file, "pz", data+5);
		RETURN(status, H5PART_SUCCESS, "H5PartReadDataFloat32");

		for (i=0; i<nparticles; i++)
		{
			FVALUE(data[6*i]  , 0.0F + (float)(i+nparticles*t), "x data");
			FVALUE(data[6*i+1], 0.1F + (float)(i+nparticles*t), "y data");
			FVALUE(data[6*i+2], 0.2F + (float)(i+nparticles*t), "z data");
			FVALUE(data[6*i+3], 0.3F + (float)(i+nparticles*t), "px data");
			FVALUE(data[6*i+4], 0.4F + (float)(i+nparticles*t), "py data");
			FVALUE(data[6*i+5], 0.5F + (float)(i+nparticles*t), "pz data");
		}

                test_read_step_attribs(file, t);
	}
}

void test_read1(void)
{
	H5PartFile *file1;

	h5part_int64_t status;

	TEST("Opening file once, read-only");
	file1 = OPEN(FILENAME,H5PART_READ);
	test_is_valid(file1);

	test_read_file_attribs(file1, 0);
	test_read_data32(file1, NPARTICLES, 1);

	status = H5PartCloseFile(file1);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
}

void test_read2(void)
{
	H5PartFile *file1;
	H5PartFile *file2;

	h5part_int64_t status;

	TEST(	"Opening file twice, read-only");
	file1 = OPEN(FILENAME,H5PART_READ);
	test_is_valid(file1);
	file2 = OPEN(FILENAME,H5PART_READ);
	test_is_valid(file2);

	test_read_strided_data32(file1, NPARTICLES, NTIMESTEPS+1);
	test_read_file_attribs(file2, 1);

	status = H5PartCloseFile(file1);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
	status = H5PartCloseFile(file2);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
}

void test_read3(void)
{
	H5PartFile *file1;

	h5part_int64_t status;

	TEST(	"Opening file once, read-only, lustre filesystem, "
		"MPI-POSIX VFD, 64KB alignment");
	file1 = OPENALIGN(FILENAME,
		H5PART_READ | H5PART_VFD_MPIPOSIX | H5PART_FS_LUSTRE,
		65536);
	test_is_valid(file1);

	TEST("Redefining step name");
	status = H5PartDefineStepName(file1, LONGNAME, 16);
	RETURN(status, H5PART_SUCCESS, "H5PartDefineStepName");

	test_read_strided_data64(file1, NPARTICLES, 0);
	test_read_file_attribs(file1, 0);

	status = H5PartCloseFile(file1);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
}

void test_read4(void)
{
	H5PartFile *file1;
	H5PartFile *file2;

	h5part_int64_t status;

	TEST(	"Opening file twice, read-only, lustre filesystem "
		"MPI-IO Independent VFD, 64K alignment");
	file1 = OPENALIGN(FILENAME,
		H5PART_READ | H5PART_VFD_MPIIO_IND | H5PART_FS_LUSTRE,
		65536);
	test_is_valid(file1);
	file2 = OPENALIGN(FILENAME,
		H5PART_READ | H5PART_VFD_MPIIO_IND | H5PART_FS_LUSTRE,
		65536);
	test_is_valid(file2);

	TEST("Redefining step name");
	status = H5PartDefineStepName(file1, LONGNAME, 16);
	RETURN(status, H5PART_SUCCESS, "H5PartDefineStepName");

	status = H5PartDefineStepName(file2, LONGNAME, 16);
	RETURN(status, H5PART_SUCCESS, "H5PartDefineStepName");

	test_read_file_attribs(file1, 1);

	status = H5PartSetStep(file2, NTIMESTEPS);
	RETURN(status, H5PART_SUCCESS, "H5PartSetStep");

	test_read_data64(file2, NPARTICLES, NTIMESTEPS-2);

	status = H5PartCloseFile(file1);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
	status = H5PartCloseFile(file2);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
}

