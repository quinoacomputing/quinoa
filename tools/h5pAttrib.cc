/* h5pAttrib.cc
   Antino Kim
   This utility will output information on h5part files accodring to the flags provided from the command line.
   The parser was imported from the example of h5dump utility with slight modifications.
*/

#include <stdio.h>
#include <stdlib.h>
#include <cctype>
#include <string.h>
#include <hdf5.h>
#include "H5Part.h"

#define MAX_LEN 100

/* Function headers */
int get_option(int argc, const char **argv, const char *opts, const struct long_options *l_opts);
static void print_help();
static void free_handler(struct arg_handler *hand, int len);
static void print_all(H5PartFile* file);
static void print_num_steps(H5PartFile* file, char * garbage);
static void print_file_attributes(H5PartFile* file, char * garbage);
static void print_step_attributes(H5PartFile* file, char *attr);
static void print_datasets(H5PartFile* file, char *attr);
static struct arg_handler* function_assign(int argc, const char *argv[]);

/* Global variables */
static int         display_all       = true;
static int         print_header      = false;
static const char* global_fname      = NULL;

/* `get_option' variables */
int         opt_err = 1;    /*get_option prints errors if this is on */
int         opt_ind = 1;    /*token pointer                          */
const char *opt_arg = NULL;        /*flag argument (or value)               */

/* indication whether the flag (option) requires an argument or not */
enum {
    no_arg = 0,         /* doesn't take an argument     */
    require_arg,        /* requires an argument	        */
};

/* struct for flags (options) */
typedef struct long_options
{
    const char  *name;          /* name of the long option              */
    int          has_arg;       /* whether we should look for an arg    */
    char         shortval;      /* the shortname equivalent of long arg
                                 * this gets returned from get_option   */
} long_options;

/* List of options in single characters */
static const char *s_opts = "hnAHa:d:";

/* List of options in full words */
static struct long_options l_opts[] =
{
    { "help", no_arg, 'h' },         // Print help page
    { "nstep", no_arg, 'n' },        // Print number of steps
    { "fileA", no_arg, 'A' },        // Print file attributes
    { "stepA", require_arg, 'a' },   // Print step attributes & values for time step n
    { "dataset", require_arg, 'd' }, // Print data sets names & values for time step n
    { "header", require_arg, 'H' },  // Print shorter version without the values
    { NULL, 0, '\0' }
};

/* a structure for handling the order command-line parameters come in */
struct arg_handler {
    void (*func)(H5PartFile *, char *);
    char *obj;
};


/************************************************************************************
***********************************  FUNCTIONS  *************************************
*************************************************************************************/


/* get_option is the parsing function that was majorly ported from h5dump utility */
int get_option(int argc, const char **argv, const char *opts, const struct long_options *l_opts)
{
    static int sp = 1;    /* character index in current token */
    int opt_opt = '?';    /* option character passed back to user */

    if (sp == 1) 
    {
        /* check for more flag-like tokens */
        if (opt_ind >= argc || argv[opt_ind][0] != '-' || argv[opt_ind][1] == '\0') 
        {
            return EOF;
        }
        else if (strcmp(argv[opt_ind], "--") == 0)
        {
            opt_ind++;
            return EOF;
        }
    }

    if (sp == 1 && argv[opt_ind][0] == '-' && argv[opt_ind][1] == '-') 
    {
        /* long command line option */
        const char *arg = &argv[opt_ind][2];
        int i;

        for (i = 0; l_opts && l_opts[i].name; i++)
        {
            size_t len = strlen(l_opts[i].name);

            if (strncmp(arg, l_opts[i].name, len) == 0)
            {
                /* we've found a matching long command line flag */
                opt_opt = l_opts[i].shortval;

                if (l_opts[i].has_arg != no_arg)
                {
                    if (arg[len] == '=')
                    {
                        opt_arg = &arg[len + 1];
                    }
                    else if (opt_ind < (argc - 1) && argv[opt_ind + 1][0] != '-')
                    {
                        opt_arg = argv[++opt_ind];
                    }
                    else if (l_opts[i].has_arg == require_arg)
                    {
                        if (opt_err)
                            fprintf(stderr, "%s: option required for \"--%s\" flag\n", argv[0], arg);

                        opt_opt = '?';
                    }
                }
                else
                {
                    if (arg[len] == '=')
                    {
                        if (opt_err)
                            fprintf(stderr, "%s: no option required for \"%s\" flag\n", argv[0], arg);

                        opt_opt = '?';
                    }

                    opt_arg = NULL;
                }

                break;
            }
        }

        if (l_opts[i].name == NULL)
        {
            /* exhausted all of the l_opts we have and still didn't match */
            if (opt_err)
                fprintf(stderr, "%s: unknown option \"%s\"\n", argv[0], arg);

            opt_opt = '?';
        }

        opt_ind++;
        sp = 1;
    }
    else
    {
        register char *cp;    /* pointer into current token */

        /* short command line option */
        opt_opt = argv[opt_ind][sp];

        if (opt_opt == ':' || (cp = strchr(opts, opt_opt)) == 0)
        {

            if (opt_err)
                fprintf(stderr, "%s: unknown option \"%c\"\n", argv[0], opt_opt);
            /* if no chars left in this token, move to next token */
            if (argv[opt_ind][++sp] == '\0')
            {
                opt_ind++;
                sp = 1;
            }

            return '?';
        }

        if (*++cp == ':')
        {

            /* if a value is expected, get it */
            if (argv[opt_ind][sp + 1] != '\0')
            {
                /* flag value is rest of current token */
                opt_arg = &argv[opt_ind++][sp + 1];
            }
            else if (++opt_ind >= argc)
            {
                if (opt_err)
                {
                    fprintf(stderr, "%s: value expected for option \"%c\"\n", argv[0], opt_opt);
                }
                opt_opt = '?';
            }
            else
            {
                /* flag value is next token */
                opt_arg = argv[opt_ind++];
            }

            sp = 1;
        }
        else 
        {
            /* set up to look at next char in token, next time */
            if (argv[opt_ind][++sp] == '\0')
            {
                /* no more in current token, so setup next token */
                opt_ind++;
                sp = 1;
            }

            opt_arg = NULL;
        }
    }

    /* return the current flag character found */
    return opt_opt;
}

/* Assigns functions according to the parsed result */
static struct arg_handler* function_assign(int argc, const char *argv[])
{
    struct arg_handler   *hand = NULL;
 
    int                  i, option;

    /* this will be plenty big enough to hold the info */
    hand = (arg_handler*)calloc((size_t)argc, sizeof(struct arg_handler));

    /* set options according to the command line */
    while ((option = get_option(argc, argv, s_opts, l_opts)) != EOF)
    {
       switch ((char)option)
       {
          case 'h': // Print help page
            print_help();
            exit(1);
          case 'A': // Print file attributes
            display_all = 0;

            for (i = 0; i < argc; i++)
            {
               if (!hand[i].func)
               {
                  hand[i].func = print_file_attributes;
                  hand[i].obj = NULL; // inserting garabage value that we won't use. (For function interface compatibility)
                  break;
               }
            }
            break;
          case 'a': // Print step attributes & values for time step n
            display_all = 0;

            for (i = 0; i < argc; i++)
            {
               if (!hand[i].func)
               {
                  hand[i].func = print_step_attributes;
                  hand[i].obj = strdup(opt_arg);
                  break;
               }
            }
            break;
          case 'd': // Print data sets names & values for time step n
            display_all = 0;

            for (i = 0; i < argc; i++)
            {
               if (!hand[i].func)
               {
                  hand[i].func = print_datasets;
                  hand[i].obj = strdup(opt_arg);
                  break;
               }
            }
            break;
          case 'n': // Print number of steps
            display_all = 0;

            for (i = 0; i < argc; i++)
            {
               if (!hand[i].func)
               {
                  hand[i].func = print_num_steps;
                  hand[i].obj = NULL; // inserting garabage value that we won't use. (For function interface compatibility)
                  break;
               }
            }
          break;
          case 'H': // Print shorter version without the values
            print_header = true;
            break;
          default:
            print_help();
            exit(1);
       }
    }
    return hand;
}

/* For printing help page */
static void
print_help () {
	fflush(stdout);
	fprintf(stdout, "\nusage: h5pAttrib [OPTIONS] file\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "  OPTIONS\n");
	fprintf(stdout, "   -h, --help           Print help page\n");
	fprintf(stdout, "   -n, --nstep          Print number of steps\n");
	fprintf(stdout, "   -A, --fileA          Print file attributes\n");
	fprintf(stdout, "   -a n, --stepA n      Print step attributes & values for time step n\n");
	fprintf(stdout, "   -d n, --dataset n    Print data sets names & values for time step n\n");
	fprintf(stdout, "   -H, --header         Print shorter version without the values\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "  Examples:\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "  1) Show file attribute names & values of sample.h5part\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "        h5pAttrib -A sample.h5part\n");
	fprintf(stdout, "\t\t\tOR\n");
	fprintf(stdout, "        h5pAttrib --fileA sample.h5part\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "  2) Show step attribute names for time step 5 of sample.h5part\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "        h5pAttrib -a 5 -H sample.h5part\n");
	fprintf(stdout, "\t\t\tOR\n");
	fprintf(stdout, "        h5pAttrib --stepA 5 -H sample.h5part\n");
	fprintf(stdout, "\n");
}


static void
print_num_steps (
	H5PartFile* file,
	char* garbage
	) {
	h5part_int64_t num_steps = H5PartGetNumSteps (file);
	fprintf(stdout, "\nNumber of steps: %lld\n", (long long)num_steps);
}


static void
print_int32 (h5part_int32_t* values, h5part_int64_t num_values) {
	h5part_int64_t k = 0;
	h5part_int64_t count = 1;
	fprintf (stdout, "\t(0): ");
	for (h5part_int64_t j = 0; j < num_values; j++) {
		fprintf (stdout, "%lld  ",(long long)values[j]);
		k++;
		if (k == 4) {
			fprintf (stdout, "\n");
			fprintf (stdout, "\t(%lld): ", (long long)count);
			k = 0;
		}
		count++;
	}
	fprintf (stdout, "\n");
}

static void
print_int64 (h5part_int64_t* values, h5part_int64_t num_values) {
	h5part_int64_t k = 0;
	h5part_int64_t count = 1;
	fprintf (stdout, "\t(0): ");
	for (h5part_int64_t j = 0; j < num_values; j++) {
		fprintf (stdout, "%lld  ",(long long)values[j]);
		k++;
		if (k == 4) {
			fprintf (stdout, "\n\t(%lld): ", (long long)count);
			k = 0;
		}
		count++;
	}
	fprintf (stdout, "\n");
}

static void
print_float64 (h5part_float64_t* values, h5part_int64_t num_values) {
	h5part_int64_t k = 0;
	h5part_int64_t count = 1;
	fprintf (stdout, "\t(0): ");
	for (h5part_int64_t j = 0; j < num_values; j++) {
		fprintf (stdout, "%lf  ", (double)values[j]);
		k++;
		if (k == 4) {
			fprintf (stdout, "\n\t(%lld): ", (long long)count);
			k = 0;
		}
		count++;
	}
	fprintf(stdout, "\n");
}

static void
print_char (char* values, h5part_int64_t num_values) {
	h5part_int64_t k = 0;
	h5part_int64_t count = 1;
	fprintf (stdout, "\t(0): ");

	for (h5part_int64_t j = 0; j < num_values; j++) {
		fprintf(stdout, "%c", values[j]);
		k++;
		if (k == 100) {
			fprintf (stdout, "\n\t(%lld): ", (long long)count);
			k = 0;
		}
	}
	fprintf(stdout, "\n");
}

static inline void
print_type (const char* s, h5part_int64_t type) {
	if (type == H5T_NATIVE_INT32) {
		fprintf(stdout, "%s type: H5T_NATIVE_INT32\n", s);
	} else if (type == H5T_NATIVE_INT64) {
		fprintf(stdout, "%s type: H5T_NATIVE_INT64\n", s);
	} else if (type == H5T_NATIVE_CHAR) {
		fprintf(stdout, "%s type: H5T_NATIVE_CHAR\n", s);
	} else if (type == H5T_NATIVE_DOUBLE) {
		fprintf(stdout, "%s type: H5T_NATIVE_DOUBLE\n", s);
	} else {
		fprintf(stdout, "%s type: UNKNOWN.\n", s);
	}
}

static inline void
print_attrib (
	H5PartFile* file,
	h5part_int64_t i,
	h5part_int64_t (*get_info)(H5PartFile*, h5part_int64_t, char*, h5part_int64_t, h5part_int64_t*, h5part_int64_t*),
	h5part_int64_t (*read_attrib)(H5PartFile*, const char*, void*)
	) {
	char* name = NULL;
	h5part_int64_t type;
	h5part_int64_t num_values;
	get_info (file, i, name, (h5part_int64_t)MAX_LEN, &type, &num_values);
	fprintf (stdout, "\tAttribute #%lld: %s\n", (long long)i, name);
	fprintf (stdout, "\tNumber of elements in attribute: %lld\n",(long long)num_values);
	print_type ("\tAttribute", type);
	if (!print_header) {
		if (type == H5T_NATIVE_INT32) {
			h5part_int32_t* values = (h5part_int32_t*)malloc(sizeof(h5part_int32_t)*num_values);
			read_attrib (file, name, values);
			print_int32 (values, num_values);
			free (values);
		} else if (type == H5T_NATIVE_INT64) {
			h5part_int64_t* values = (h5part_int64_t*)malloc(sizeof(h5part_int64_t)*num_values);
			read_attrib (file, name, values);
			print_int64 (values, num_values);
			free (values);
		} else if (type == H5T_NATIVE_CHAR) {
			char* values = (char*)malloc(sizeof(char)*num_values);
			read_attrib (file, name, values);
			print_char (values, num_values);
			free (values);
		} else if (type == H5T_NATIVE_DOUBLE) {
			h5part_float64_t* values = (h5part_float64_t*)malloc(sizeof(h5part_float64_t)*num_values);
			read_attrib (file, name, values);
			print_float64 (values, num_values);
			free (values);
		}
		fprintf(stdout, "\n");
	}
}

static void
print_file_attribute (
	H5PartFile* file,
	h5part_int64_t i
	) {
	print_attrib (file, i, H5PartGetFileAttribInfo, H5PartReadFileAttrib);
}

static void
print_file_attributes (
	H5PartFile* file,
	char * garbage
	) {
	h5part_int64_t num_attribs = H5PartGetNumFileAttribs(file);
	fprintf(stdout, "Number of file attributes: %lld\n", (long long)num_attribs);

	h5part_int64_t i;
	for (i = 0; i < num_attribs; i++) {
		print_file_attribute (file, i);
	}
}

static void
print_step_attribute (
	H5PartFile* file,
	h5part_int64_t i
	) {
	print_attrib (file, i, H5PartGetStepAttribInfo, H5PartReadStepAttrib);
}

/* For printing dataset names [& values] */
static void
_print_step_attributes (
	H5PartFile* file,
	h5part_int64_t step
	) {
	h5part_int64_t num_attribs = H5PartGetNumStepAttribs(file);
	fprintf(stdout, "Number of step attributes in step #%lld: %lld\n", step, (long long)num_attribs);

	h5part_int64_t i;
	for (i = 0; i < num_attribs; i++) {
		print_step_attribute (file, i);
	}
}


static void
print_step_attributes(
	H5PartFile* file,
	char* step_sz
	) {
	h5part_int64_t step = atoi(step_sz);
	H5PartSetStep (file, step);
	_print_step_attributes (file, step);
}

static inline void
print_dataset (
	H5PartFile* file, h5part_int64_t i
	) {
	char name[MAX_LEN];
	h5part_int64_t type;
	h5part_int64_t num_values;

	H5PartGetDatasetInfo (file, i, name, MAX_LEN, &type, &num_values);
	fprintf (stdout, "\n\tDataset name #%lld: %s\n",(long long)i, name);
	print_type ("\tDataset", type);
	fprintf (stdout, "\tNumber of elements in the dataset: %lld\n", (long long)num_values);
	if (!print_header) {
		if (type == H5T_NATIVE_INT32) {
			h5part_int32_t* values = (h5part_int32_t*)malloc(sizeof(*values)*num_values);
			H5PartReadDataInt32 (file, name, values);
			print_int32 (values, num_values);
			free (values);
		} else if (type == H5T_NATIVE_INT64) {
			h5part_int64_t* values = (h5part_int64_t*)malloc(sizeof(*values)*num_values);
			H5PartReadDataInt64 (file, name, values);
			print_int64 (values, num_values);
			free (values);
		} else if (type == H5T_NATIVE_DOUBLE) {
			h5part_float64_t* values = (h5part_float64_t*)malloc (sizeof(*values)*num_values);
			H5PartReadDataFloat64 (file, name, values);
			print_float64 (values, num_values);
			free (values);
		} else {
			fprintf(stdout, "Dataset cannot be read.\n");
		}
		fprintf(stdout, "\n");
         }
}

/* For printing dataset names [& values] */
static void
_print_datasets (
	H5PartFile* file,
	h5part_int64_t step
	) {
	h5part_int64_t num_datasets = H5PartGetNumDatasets (file);

	fprintf (stdout, "Number of datasets in step #%lld: %lld\n", step, num_datasets);

	h5part_int64_t i;
	for (i = 0; i < num_datasets; i++) {
		print_dataset (file, i);
	}
}

static void
print_datasets (
	H5PartFile* file,
	char* step_sz
	) {
	h5part_int64_t step = atoi(step_sz);
	H5PartSetStep (file, step);
	_print_datasets (file, step);
}

/* For priting everything (default option when no flags provided.) */
static void
print_all(
	H5PartFile* file
	) {
	fprintf(stdout, "\nDump result for \"%s\"...\n\n", global_fname);
	print_num_steps (file, NULL);
	print_file_attributes (file, NULL);
	int num_steps = H5PartGetNumSteps(file);
	for (h5part_int64_t step = 0; step < num_steps; step++) {
		H5PartSetStep (file, step);
		_print_step_attributes (file, step);
		_print_datasets (file, step);
	}
}


/* Frees argument handlers */
static void
free_handler (
	struct arg_handler *hand,
	int len
	) {
	int i;

	for (i = 0; i < len; i++) {
		free(hand[i].obj);
	}
	
	free(hand);
}

int
main (
	int argc,
	const char *argv[]
	) {
	/* Numerous variables */
	struct arg_handler* hand = NULL;
	int i;
	H5PartFile* h5file = NULL;
	const char* fname = NULL;

	//h5pAttrib_function_table = &major_function_table;
	/* Take care of the command line options */
	hand = function_assign (argc, argv);
	
	if (argc <= opt_ind) {
		fprintf (stdout, "missing file name\n");
		print_help ();
		exit(1);
	}

	/* Check for conflicting options */
	/* There are none. If -H is appended to non-compatible flags, just ignore. */

	/* Process accordingly */
	fname = argv[opt_ind];
	global_fname = fname; // To use in funtions
	h5file = H5PartOpenFile (fname,H5PART_READ);

	if (h5file == NULL) {
		fprintf(stdout, "unable to open file %s\n", fname);
		print_help ();
		exit(1);
	}

	if (display_all) {
		print_all (h5file);
	} else {
		for (i = 0; i < argc; i++) {
			if (hand[i].func) {
				hand[i].func (h5file, hand[i].obj);
			}
		}
	}

	free_handler (hand, argc);
	H5PartCloseFile (h5file);
	fprintf (stdout,"done\n");
	return 0;
}
