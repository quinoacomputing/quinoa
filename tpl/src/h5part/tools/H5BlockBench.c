/* 
 * File:        H5BlockBench.c
 * Author:      Mark Howison
 * Created:     10-17-2008
 * Description: Benchmark application for H5Block, with similar
 *              functionality to IOR.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <H5Part.h>
#include <H5Block.h>

#define VERBOSITY_LOW     1
#define VERBOSITY_MEDIUM  3
#define VERBOSITY_HIGH    5

#define ONE_MEGABYTE 1048576.0

#ifdef PARALLEL_IO

/* Parameters ***********************************************************/
struct Params {
int         rank;               // MPI
int         procs;              // MPI
int         segments;           // number of segments (i.e. time steps)
int         repetitions;        // number of times to repeat the test
int         read;               // enable read test
int         write;              // enable write test
int         nofillvals;         // disable fill values in HDF5
int         verbosity;          // verbosity level
char*       filename;           // path to the test file
char        flags;              // file open flags
double      aggregate_size;     // aggregate size in MB
h5part_int64_t alignment;       // HDF5 alignment in bytes
h5part_int64_t block_size;      // size of the block data
h5part_int64_t bdims[3];        // block dimensions
h5part_int64_t cdims[3];        // chunk dimensions
h5part_int64_t layout[6];       // H5Block layout
};
typedef struct Params Params;


/************************************************************************/
/* Argument helper functions                                            */
/************************************************************************/

void print_usage ()
{
    printf ("Usage:\n");
    printf (" -a\talignment in bytes [1048576]\n");
    printf (" -b\tblock dimensions [64 64 64]\n");
    printf (" -c\tchunk dimentions [0 0 0]\n");
    printf (" -i\trepetitions [1]\n");
    printf (" -n\tHDF5 no fill [0]\n");
    printf (" -o\toutput dir [output]\n");
    printf (" -r\tperform read test [0]\n");
    printf (" -s\tsegments [1]\n");
    printf (" -v\tverbosity level [1]\n");
    printf (" -w\tpreform write test [0]\n");
    printf (" -lustre\tenable lustre-specific tuning\n");
    printf (" -posix\tuse the MPI-POSIX VFD\n");
}

void parse_args (int argc, char** argv, Params* p)
{
    int i;
    
    if (argc < 2) {
        if (p->rank == 0) print_usage();
        MPI_Finalize();
        exit (EXIT_SUCCESS);
    }

    // default values
    p->bdims[0]    = 64;
    p->bdims[1]    = 64;
    p->bdims[2]    = 64;
    p->cdims[0]    = 0;
    p->cdims[1]    = 0;
    p->cdims[2]    = 0;
    p->segments    = 1;
    p->repetitions = 1;
    p->read        = 0;
    p->write       = 0;
    p->alignment   = 1048576;
    p->nofillvals  = 0;
    p->verbosity   = VERBOSITY_MEDIUM;
    p->filename    = "output";
    p->flags       = H5PART_WRITE;
   
    i = 1;
    while (i < argc)
    {
        // block dimensions
        if (strcmp(argv[i],"-b") == 0)
        {
            i++;
            p->bdims[0] = atoi (argv[i]);
            i++;
            p->bdims[1] = atoi (argv[i]);
            i++;
            p->bdims[2] = atoi (argv[i]);
        }
        // chunk dimensions
        else if (strcmp(argv[i],"-c") == 0)
        {
            i++;
            p->cdims[0] = atoi (argv[i]);
            i++;
            p->cdims[1] = atoi (argv[i]);
            i++;
            p->cdims[2] = atoi (argv[i]);
        }
        // segments
        else if (strcmp(argv[i],"-s") == 0)
        {
            i++;
            p->segments = atoi (argv[i]);
        }
        // repetitions
        else if (strcmp(argv[i],"-i") == 0)
        {
            i++;
            p->repetitions = atoi (argv[i]);
        }
        // read
        else if (strcmp(argv[i],"-r") == 0)
        {
            p->read = 1;
        }
        // write
        else if (strcmp(argv[i],"-w") == 0)
        {
            p->write = 1;
        }
        // alignment
        else if (strcmp(argv[i],"-a") == 0)
        {
            i++;
            p->alignment = atoi (argv[i]);
        }
        // nofillvals 
        else if (strcmp(argv[i],"-n") == 0)
        {
            p->nofillvals = 1;
        }
        // verbosity level
        else if (strcmp(argv[i],"-v") == 0)
        {
            i++;
            p->verbosity = atoi (argv[i]);
        }
        // filename
        else if (strcmp(argv[i],"-o") == 0)
        {
            i++;
            p->filename = (char*) malloc (strlen (argv[i]) + 1);
            strcpy (p->filename, argv[i]);
        }
        else if (strcmp(argv[i],"-lustre") == 0)
        {
            p->flags |= H5PART_FS_LUSTRE;
        }
        else if (strcmp(argv[i],"-posix") == 0)
        {
            p->flags |= H5PART_VFD_MPIPOSIX;
        }
        // print usage
        else if (strcmp(argv[i],"--help") == 0)
        {
            if (p->rank == 0) print_usage();
            MPI_Finalize();
            exit (EXIT_SUCCESS);
        }
        else
        {
            if (p->rank == 0) {
                fprintf (stderr, "%s: unrecognized argument %s \n",
                        argv[0], argv[i]);
                print_usage();
            }
            MPI_Finalize();
            exit (EXIT_FAILURE);
        }
    i++;
    }
}


/************************************************************************/
/* Layout functions                                                     */
/************************************************************************/

int nth_root_int_divisor (const int m, const int n)
{
    int i, root;
    double p;

    p = 1.0 / (double) n;
    root = (int) ceil ( pow ((double) m, p) );
    for (i=root; i<=m; i++)
    {
        if (m % i == 0) return i;
    }

    return i;
}

void set_layout (Params* p)
{
   int i, j, k;
   int x, y, z;

   x = nth_root_int_divisor (p->procs, 3);
   y = nth_root_int_divisor (p->procs / x, 2);
   z = p->procs / x / y;

   if (p->verbosity >= VERBOSITY_HIGH && p->rank == 0) {
       printf ("creating (%d,%d,%d) layout\n", x, y, z);
   }

   i = p->rank % x;
   j = (p->rank / x) % y;
   k = p->rank / (x * y);

   p->layout[0] =    i    * p->bdims[0];
   p->layout[1] = (i + 1) * p->bdims[0] - 1;
   p->layout[2] =    j    * p->bdims[1];
   p->layout[3] = (j + 1) * p->bdims[1] - 1;
   p->layout[4] =    k    * p->bdims[2];
   p->layout[5] = (k + 1) * p->bdims[2] - 1;

   if (p->verbosity >= VERBOSITY_HIGH) {
       printf ("rank %d: (%d,%d,%d) [%lld,%lld]x[%lld,%lld]x[%lld,%lld]\n",
               p->rank, i, j, k,
               (long long)p->layout[0], (long long)p->layout[1],
               (long long)p->layout[2], (long long)p->layout[3],
               (long long)p->layout[4], (long long)p->layout[5]);
   }
}

void check_cdims (Params* p)
{
    int i;

    for (i=0; i<3; i++) {
        if (p->bdims[i] % p->cdims[i] != 0) {
            if (p->rank == 0) {
                fprintf (stderr,
                        "Chunk dim %d does not divide block dim %d!\n",
                        i, i);
            }
            MPI_Barrier (MPI_COMM_WORLD);
            exit (EXIT_FAILURE);
        }
    }
}


/************************************************************************/
/* Data functions                                                       */
/************************************************************************/

void create_data (float* data, Params* p)
{
    int i;

    if (p->verbosity >= VERBOSITY_HIGH) {
        printf ("rank %d: creating random data\n", p->rank);
    }

    for (i=0; i<p->block_size; i++) {
        data[i] = (float) random();
    }
}

double write_data (float* data, int iter, Params* p)
{
    int             i;
    double          start_time;
    double          open_time;
    double          write_time;
    double          close_time;
    double          total_time;
    double          sum_time;
    double          open_mean;
    double          write_mean;
    double          close_mean;
    double          bandwidth;
    float*          segment;
    char*           filename;
    H5PartFile*     file;
    h5part_int64_t  status;
    
    if (p->verbosity >= VERBOSITY_HIGH) {
        printf ("rank %d: writing data\n", p->rank);
    }

    start_time = MPI_Wtime();
    
    filename = (char*) malloc (strlen (p->filename) + 64);
    sprintf (filename, "%s/%d.h5", p->filename, iter);

    file = H5PartOpenFileParallelAlign (filename,
            p->flags, MPI_COMM_WORLD, p->alignment);
    if (!file) {
        fprintf (stderr,
                "rank %d: could not open H5Part file!\n", p->rank);
        MPI_Barrier (MPI_COMM_WORLD);
        exit (EXIT_FAILURE);
    }

    open_time = MPI_Wtime() - start_time;

    if (p->cdims[0] > 0 &&  p->cdims[1] > 0 && p->cdims[2] > 0) {
        status = H5BlockDefine3DChunkDims (file,
                p->cdims[0], p->cdims[1], p->cdims[2]);
        if (status != H5PART_SUCCESS) {
            fprintf (stderr,
                    "rank %d: H5Block chunk error!", p->rank);
        }
    }

    status = H5BlockDefine3DFieldLayout (file,
            p->layout[0], p->layout[1],
            p->layout[2], p->layout[3],
            p->layout[4], p->layout[5]);
    if (status != H5PART_SUCCESS) {
        fprintf (stderr,
                "rank %d: H5Block layout error!", p->rank);
    }

    segment = data;

    for (i=0; i<p->segments; i++) {
        status = H5PartSetStep (file, i);
        if (status != H5PART_SUCCESS) {
            fprintf (stderr, "rank %d: H5PartSetStep error!", p->rank);
        }

        status = H5Block3dWriteScalarFieldFloat32 (file, "test", segment);
        if (status != H5PART_SUCCESS) {
            fprintf (stderr, "rank %d: H5Block write error!", p->rank);
        }

        segment += p->block_size;
    }

    write_time = (MPI_Wtime() - start_time) - open_time;
    
    H5PartCloseFile (file);
    
    close_time = (MPI_Wtime() - start_time) - write_time - open_time;

    total_time = open_time + write_time + close_time;

    if (p->verbosity >= VERBOSITY_HIGH) {
        printf ("rank %d: write\t%.3f\t%.3f\t%.3f\t%.3f\n", p->rank,
                open_time, write_time, close_time, total_time);
    }

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Reduce (&open_time, &sum_time, 1, MPI_DOUBLE,
            MPI_SUM, 0, MPI_COMM_WORLD);
    open_mean = sum_time / p->procs;
    MPI_Reduce (&write_time, &sum_time, 1, MPI_DOUBLE,
            MPI_SUM, 0, MPI_COMM_WORLD);
    write_mean = sum_time / p->procs;
    MPI_Reduce (&close_time, &sum_time, 1, MPI_DOUBLE,
            MPI_SUM, 0, MPI_COMM_WORLD);
    close_mean = sum_time / p->procs;
    bandwidth = p->aggregate_size / total_time;

    if (p->verbosity >= VERBOSITY_MEDIUM && p->rank == 0) {
        printf ("write\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\n", bandwidth,
                open_mean, write_mean, close_mean, total_time);
    }

    return bandwidth;
}

 
/************************************************************************/
/* Main procedure                                                       */
/************************************************************************/

int main (int argc, char** argv)
{
    int         i;
    Params      p;
    float*      data;
    double      chunk_size;
    double      data_size;
    double      bandwidth;
    double      write_max;
    double      write_min;
    double      write_mean;
    double      read_max;
    double      read_min;
    double      read_mean;
    time_t      rawtime;
    struct tm * timeinfo;

    // initialize MPI
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &p.rank);
    MPI_Comm_size (MPI_COMM_WORLD, &p.procs);
    
    parse_args (argc, argv, &p);
    //check_cdims (&p);

    H5PartSetVerbosityLevel (p.verbosity);

    if (p.verbosity >= VERBOSITY_MEDIUM && p.rank == 0) {
        time (&rawtime);
        timeinfo = localtime (&rawtime);
        printf ("Started: %s", asctime (timeinfo));
        printf ("Command line:\n");
        for (i=0; i<argc; i++) {
            printf ("%s ",argv[i]);
        }
        printf ("\n");
    }

    p.block_size = p.bdims[0]*p.bdims[1]*p.bdims[2];
    chunk_size = p.cdims[0]*p.cdims[1]*p.cdims[2] * sizeof(float);
    chunk_size /= ONE_MEGABYTE;
    data_size = p.segments * p.block_size * sizeof(float);
    data_size /= ONE_MEGABYTE;
    p.aggregate_size = p.procs * data_size;

    set_layout (&p);

    if (p.verbosity >= VERBOSITY_MEDIUM && p.rank == 0) {
        printf ("chunk_size\t\t= %.1f MB\n", chunk_size);
        printf ("block_size\t\t= %.1f MB\n",
                p.block_size * sizeof(float) / 1048576.0);
        printf (" x %5d segment(s)\t= %.1f MB per proc\n",
                p.segments, data_size);
        printf (" x %5d procs\t\t= %.1f MB aggregate\n",
                p.procs, p.aggregate_size);
    }

    if (p.verbosity >= VERBOSITY_HIGH) {
        printf ("rank %d: mallocing block data\n", p.rank);
    }

    data = (float*) malloc (p.segments * p.block_size * sizeof(float));
    if (!data) {
        fprintf (stderr,
                "rank %d: could not malloc block data!\n", p.rank);
        MPI_Barrier (MPI_COMM_WORLD);
        exit (EXIT_FAILURE);
    }

    if (p.verbosity >= VERBOSITY_MEDIUM && p.rank == 0) {
        printf ("access\tbw\topen\twr/rd\tclose\ttotal\n");
        printf ("(wr/rd)\t(MB/s)\t(s)\t(s)\t(s)\t(s)\n");
        printf ("------\t--\t----\t-----\t-----\t-----\n");
    }

    write_max = 0.0;
    write_min = 10e300;
    for (i=0; i<p.repetitions; i++) {
        create_data (data, &p);
        MPI_Barrier (MPI_COMM_WORLD);
        if (p.write) {
            bandwidth = write_data (data, i, &p);
            if (bandwidth > write_max) write_max = bandwidth;
            if (bandwidth < write_min) write_min = bandwidth;
            write_mean += bandwidth;
        }
    }
    write_mean /= p.repetitions;

   if (p.verbosity >= VERBOSITY_LOW && p.rank == 0) {
        printf ("********************************\n");
        printf ("** Aggregate Bandwidth (MB/s) **\n");
        printf ("access  max     min     mean    \n");
        printf ("------  ---     ---     ----    \n");
        if (p.write) {
            printf ("write\t%.1f\t%.1f\t%.1f\n",
                    write_max, write_min, write_mean);
        }
        if (p.read) {
            printf ("read\t%.1f\t%.1f\t%.1f\n",
                    read_max, read_min, read_mean);
        }
        printf ("********************************\n");
        time (&rawtime);
        timeinfo = localtime (&rawtime);
        printf ("Finished: %s", asctime (timeinfo));
    }
    
    if (p.verbosity >= VERBOSITY_HIGH) {
        printf ("rank %d: Exiting!\n", p.rank);
    }
            
    MPI_Finalize();
    return (EXIT_SUCCESS);
}

#else
#error This file only works when PARALLEL_IO is enabled.
#endif
