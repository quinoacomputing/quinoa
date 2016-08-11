/* 
 * File:        H5PartBench.c
 * Author:      Mark Howison
 * Created:     11-29-2008
 * Description: Benchmark application for H5Part, with similar
 *              functionality to IOR.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <H5Part.h>

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
h5part_int64_t particles;       // number of particles to write per block
h5part_int64_t blocks;          // number of blocks to write per segment
};
typedef struct Params Params;


/************************************************************************/
/* Argument helper functions                                            */
/************************************************************************/

void print_usage ()
{
    printf ("Usage:\n");
    printf (" -a\talignment in bytes\n");
    printf (" -b\tblocks per segment (i.e. chunks per timestep)\n");
    printf (" -i\trepetitions\n");
    printf (" -n\tHDF5 no fill\n");
    printf (" -o\toutput dir\n");
    printf (" -p\tparticles per block\n");
    printf (" -r\tperform read test\n");
    printf (" -s\tsegments (i.e. H5Part timesteps)\n");
    printf (" -v\tverbosity level\n");
    printf (" -w\tpreform write test\n");
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
    p->blocks      = 1;
    p->particles   = 262144;
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
        // blocks
        if (strcmp(argv[i],"-b") == 0)
        {
            i++;
            p->blocks = atoi (argv[i]);
        }
        // particles
        else if (strcmp(argv[i],"-p") == 0)
        {
            i++;
            p->particles = atoi (argv[i]);
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
/* Data functions                                                       */
/************************************************************************/

void create_data (float* data, Params* p)
{
    int i;

    if (p->verbosity >= VERBOSITY_HIGH) {
        printf ("rank %d: creating random data\n", p->rank);
    }

    for (i=0; i<(p->blocks * p->particles); i++) {
        data[i] = (float) random();
    }
}

double write_data (float* data, int iter, Params* p)
{
    int             i,j;
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
    char            var_name[64];
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

    H5PartSetNumParticles (file, p->particles);
    if (p->verbosity >= VERBOSITY_HIGH) {
        printf ("rank %d: %ld particles\n", p->rank, p->particles);
    }

    open_time = MPI_Wtime() - start_time;

    segment = data;

    for (i=1; i<=p->segments; i++) {
        status = H5PartSetStep (file, i);
        if (status != H5PART_SUCCESS) {
            fprintf (stderr, "rank %d: H5PartSetStep error!", p->rank);
        }
        for (j=0; j<p->blocks; j++) {
            sprintf (var_name, "test%d", j);
            status = H5PartWriteDataFloat32 (file, var_name, segment);
            if (status != H5PART_SUCCESS) {
                fprintf (stderr,
                        "rank %d: H5PartWriteDataFloat32 error!", p->rank);
            }
            segment += p->particles;
        }
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

    free (filename);

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
    double      block_size;
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

    block_size =  p.particles * sizeof(float);
    block_size /= ONE_MEGABYTE;
    data_size = p.segments * p.blocks * block_size;
    p.aggregate_size = p.procs * data_size;

    if (p.verbosity >= VERBOSITY_MEDIUM && p.rank == 0) {
        printf ("block_size\t\t= %.1f MB\n", block_size);
        printf ("segment_size\t\t= %.1f MB\n",
                p.blocks * block_size);
        printf (" x %5d segment(s)\t= %.1f MB per proc\n",
                p.segments, data_size);
        printf (" x %5d procs\t\t= %.1f MB aggregate\n",
                p.procs, p.aggregate_size);
    }

    if (p.verbosity >= VERBOSITY_HIGH) {
        printf ("rank %d: mallocing block data\n", p.rank);
    }

    data = (float*) malloc (p.segments * p.blocks * p.particles
            * sizeof(float));
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
