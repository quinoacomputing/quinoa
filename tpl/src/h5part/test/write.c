#include <H5Part.h>

#include "testframe.h"
#include "params.h"

static void
test_write_file_attribs(H5PartFile *file, int position)
{
	h5part_int64_t status;
	char name[ATTR_NAME_SIZE];

	TEST("Writing file attributes");

	get_attr_name(name, "str", position);
	status = H5PartWriteFileAttribString(file, name, ATTR_STR_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteFileAttribString");

	get_attr_name(name, "i32", position);
	status = H5PartWriteFileAttribInt32(file, name, ATTR_INT32_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteFileAttribInt32");

	get_attr_name(name, "i64", position);
	status = H5PartWriteFileAttribInt64(file, name, ATTR_INT64_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteFileAttribInt64");

	get_attr_name(name, "f32", position);
	status = H5PartWriteFileAttribFloat32(file, name, ATTR_FLOAT_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteFileAttribFloat32");

	get_attr_name(name, "f64", position);
	status = H5PartWriteFileAttribFloat64(file, name, ATTR_FLOAT_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteFileAttribFloat64");
}

static void
test_write_step_attribs(H5PartFile *file, int position)
{
	h5part_int64_t status;
	char name[ATTR_NAME_SIZE];

	TEST("Writing step attributes");

	get_attr_name(name, "str", position);
	status = H5PartWriteStepAttribString(file, name, ATTR_STR_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteStepAttribString");

	get_attr_name(name, "i32", position);
	status = H5PartWriteStepAttribInt32(file, name, ATTR_INT32_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteStepAttribInt32");

	get_attr_name(name, "i64", position);
	status = H5PartWriteStepAttribInt64(file, name, ATTR_INT64_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteStepAttribInt64");

	get_attr_name(name, "f32", position);
	status = H5PartWriteStepAttribFloat32(file, name, ATTR_FLOAT_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteStepAttribFloat32");

	get_attr_name(name, "f64", position);
	status = H5PartWriteStepAttribFloat64(file, name, ATTR_FLOAT_VAL);
	RETURN(status, H5PART_SUCCESS, "H5PartWriteStepAttribFloat64");
}

static void
test_write_data64(H5PartFile *file, int nparticles, int step)
{
	int i,t;
	h5part_int64_t status, val;

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

	/* invalid stride will produce a warning */
	status = H5PartSetNumParticlesStrided(file, nparticles, -1);
	RETURN(status, H5PART_SUCCESS, "H5PartSetNumParticlesStrided");

	/* invalid nparticles will produce an error */
	status = H5PartSetNumParticlesStrided(file, -1, 2);
	RETURN(status, H5PART_ERR_INVAL, "H5PartSetNumParticlesStrided");

#if PARALLEL_IO
	TEST("Setting throttle");
	status = H5PartSetThrottle(file, 2);
	RETURN(status, H5PART_SUCCESS, "H5PartSetThrottle");
#endif

	TEST("Writing 64-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (i=0; i<nparticles; i++)
		{
			x[i]  = 0.0 + (double)(i+nparticles*t);
			y[i]  = 0.1 + (double)(i+nparticles*t);
			z[i]  = 0.2 + (double)(i+nparticles*t);
			px[i] = 0.3 + (double)(i+nparticles*t);
			py[i] = 0.4 + (double)(i+nparticles*t);
			pz[i] = 0.5 + (double)(i+nparticles*t);
			id[i] = i + nparticles*t;
		}

		val = H5PartHasStep(file, t);

		status = H5PartSetStep(file, t);
		RETURN(status, H5PART_SUCCESS, "H5PartSetStep");

		if (val == 0) test_write_step_attribs(file, t);

		status = H5PartSetNumParticles(file, nparticles);
		RETURN(status, H5PART_SUCCESS, "H5PartSetNumParticles");

		status = H5PartWriteDataFloat64(file, "x", x);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "y", y);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "z", z);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");
		
		status = H5PartWriteDataFloat64(file, "px", px);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "py", py);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "pz", pz);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");
		
		status = H5PartWriteDataInt64(file, "id", id);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataInt64");
	}
}

static void
test_write_strided_data64(H5PartFile *file, int nparticles, int step)
{
	int i,t;
	h5part_int64_t status;

	double *data;

	data=(double*)malloc(6*nparticles*sizeof(double));

	status = H5PartSetNumParticlesStrided(file, nparticles, 6);
	RETURN(status, H5PART_SUCCESS, "H5PartSetNumParticlesStrided");

	TEST("Writing 64-bit strided data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (i=0; i<nparticles; i++)
		{
			data[6*i]   = 0.0 + (double)(i+nparticles*t);
			data[6*i+1] = 0.1 + (double)(i+nparticles*t);
			data[6*i+2] = 0.2 + (double)(i+nparticles*t);
			data[6*i+3] = 0.3 + (double)(i+nparticles*t);
			data[6*i+4] = 0.4 + (double)(i+nparticles*t);
			data[6*i+5] = 0.5 + (double)(i+nparticles*t);
		}

		status = H5PartSetStep(file, t);
		RETURN(status, H5PART_SUCCESS, "H5PartSetStep");

		status = H5PartSetNumParticlesStrided(file, nparticles, 6);
		RETURN(status, H5PART_SUCCESS, "H5PartSetNumParticlesStrided");

		status = H5PartWriteDataFloat64(file, "x", data+0);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "y", data+1);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "z", data+2);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");
		
		status = H5PartWriteDataFloat64(file, "px", data+3);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "py", data+4);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");

		status = H5PartWriteDataFloat64(file, "pz", data+5);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat64");

                test_write_step_attribs(file, t);
	}
}

static void
test_write_data32(H5PartFile *file, int nparticles, int step)
{
	int i,t;
	h5part_int32_t status, val;
	int rank, nprocs;

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

	status = H5PartSetNumParticles(file, nparticles);
	RETURN(status, H5PART_SUCCESS, "H5PartSetNumParticles");

#if PARALLEL_IO
	/* will generate a warning since we are in MPI-IO Collective mode */
	TEST("Setting throttle");
	status = H5PartSetThrottle(file, 2);
	RETURN(status, H5PART_SUCCESS, "H5PartSetThrottle");

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#else
	rank = 0;
	nprocs = 1;
#endif

	TEST("Writing 32-bit data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (i=0; i<nparticles; i++)
		{
			x[i]  = 0.0F + (float)(i+nparticles*t);
			y[i]  = 0.1F + (float)(i+nparticles*t);
			z[i]  = 0.2F + (float)(i+nparticles*t);
			px[i] = 0.3F + (float)(i+nparticles*t);
			py[i] = 0.4F + (float)(i+nparticles*t);
			pz[i] = 0.5F + (float)(i+nparticles*t);
			id[i] = i + nparticles*t;
		}

		val = H5PartHasStep(file, t);
		if (val == 0) {
			status = H5PartSetStep(file, t);
			RETURN(status, H5PART_SUCCESS, "H5PartSetStep");
		}

		/* test a two-part write using views */
		status = H5PartSetView(file,
			rank*nparticles,
			rank*nparticles + 31);
		RETURN(status, H5PART_SUCCESS, "H5PartSetView");

		status = H5PartWriteDataFloat32(file, "x", x);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

                test_write_step_attribs(file, t);

		status = H5PartWriteDataFloat32(file, "y", y);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "z", z);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");
		
		status = H5PartWriteDataFloat32(file, "px", px);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "py", py);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "pz", pz);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");
		
		status = H5PartWriteDataInt32(file, LONGNAME, id);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataInt32");

		/* the second write phase... */
		status = H5PartSetView(file,
			rank*nparticles + 32,
			rank*nparticles + nparticles - 1);
		RETURN(status, H5PART_SUCCESS, "H5PartSetView");
		/* offset the input arrays */
		status = H5PartWriteDataFloat32(file, "x", x+32);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "y", y+32);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "z", z+32);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");
		
		status = H5PartWriteDataFloat32(file, "px", px+32);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "py", py+32);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "pz", pz+32);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");
		
		status = H5PartWriteDataInt32(file, LONGNAME, id+32);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataInt32");
	}
}

static void
test_write_strided_data32(H5PartFile *file, int nparticles, int step)
{
	int i,t;
	h5part_int64_t status;

	float *data;

	data=(float*)malloc(6*nparticles*sizeof(float));

	TEST("Writing 32-bit strided data");

	for (t=step; t<step+NTIMESTEPS; t++)
	{
		for (i=0; i<nparticles; i++)
		{
			data[6*i]   = 0.0F + (float)(i+nparticles*t);
			data[6*i+1] = 0.1F + (float)(i+nparticles*t);
			data[6*i+2] = 0.2F + (float)(i+nparticles*t);
			data[6*i+3] = 0.3F + (float)(i+nparticles*t);
			data[6*i+4] = 0.4F + (float)(i+nparticles*t);
			data[6*i+5] = 0.5F + (float)(i+nparticles*t);
		}

		status = H5PartSetStep(file, t);
		RETURN(status, H5PART_SUCCESS, "H5PartSetStep");

                test_write_step_attribs(file, t);

		status = H5PartSetNumParticlesStrided(file, nparticles, 6);
		RETURN(status, H5PART_SUCCESS, "H5PartSetNumParticlesStrided");

		status = H5PartWriteDataFloat32(file, "x", data+0);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "y", data+1);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "z", data+2);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");
		
		status = H5PartWriteDataFloat32(file, "px", data+3);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "py", data+4);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");

		status = H5PartWriteDataFloat32(file, "pz", data+5);
		RETURN(status, H5PART_SUCCESS, "H5PartWriteDataFloat32");
	}
}

void test_write1(void)
{
	H5PartFile *file1;

	h5part_int64_t status;

	TEST("Opening file once, write-truncate");
	file1 = OPEN(FILENAME,H5PART_WRITE);
	test_is_valid(file1);

	test_write_data32(file1, NPARTICLES, 1);
	test_write_file_attribs(file1, 0);

	status = H5PartCloseFile(file1);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
}

void test_write2(void)
{
	H5PartFile *file1;
	H5PartFile *file2;

	h5part_int64_t status;

	TEST("Opening file twice, write-append + read-only");
	file1 = OPEN(FILENAME,H5PART_APPEND);
	test_is_valid(file1);
	file2 = OPEN(FILENAME,H5PART_READ);
	test_is_valid(file2);

	test_write_strided_data32(file1, NPARTICLES, NTIMESTEPS+1);
	test_write_file_attribs(file1, 1);

	status = H5PartCloseFile(file1);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
	status = H5PartCloseFile(file2);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
}

void test_write3(void)
{
	H5PartFile *file1;

	h5part_int64_t status;

	TEST(	"Opening file once, write-truncate, lustre filesyste, "
		"MPI-POSIX VFD, 1KB alignment");
	file1 = OPENALIGN(FILENAME,
		H5PART_WRITE | H5PART_VFD_MPIPOSIX | H5PART_FS_LUSTRE,
		1024);
	test_is_valid(file1);

	TEST("Redefining step name");
	status = H5PartDefineStepName(file1, LONGNAME, 16);
	RETURN(status, H5PART_SUCCESS, "H5PartDefineStepName");

	test_write_strided_data64(file1, NPARTICLES, 0);
	test_write_file_attribs(file1, 0);

	status = H5PartCloseFile(file1);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
}

void test_write4(void)
{
	H5PartFile *file1;
	H5PartFile *file2;

	h5part_int64_t status;

	TEST(	"Opening file twice, write-append + read-only, "
		"lustre filesystem, MPI-IO Independent VFD, 1K alignment");
	file1 = OPENALIGN(FILENAME,
		H5PART_APPEND | H5PART_VFD_MPIIO_IND | H5PART_FS_LUSTRE,
		1024);
	test_is_valid(file1);
	file2 = OPENALIGN(FILENAME,
		H5PART_READ | H5PART_VFD_MPIIO_IND | H5PART_FS_LUSTRE,
		1024);
	test_is_valid(file2);

	TEST("Redefining step name");
	status = H5PartDefineStepName(file1, LONGNAME, 16);
	RETURN(status, H5PART_SUCCESS, "H5PartDefineStepName");

	status = H5PartDefineStepName(file2, LONGNAME, 16);
	RETURN(status, H5PART_SUCCESS, "H5PartDefineStepName");

	status = H5PartSetChunkSize(file1, NPARTICLES);
	RETURN(status, H5PART_SUCCESS, "H5PartSetChunkSize");

	test_write_data64(file1, NPARTICLES, NTIMESTEPS-2);
	test_write_file_attribs(file1, 1);

	status = H5PartCloseFile(file1);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
	status = H5PartCloseFile(file2);
	RETURN(status, H5PART_SUCCESS, "H5PartCloseFile");
}

