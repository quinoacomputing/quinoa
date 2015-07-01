 
/* ufile.h for ANSI C */
#ifndef UFILE_H
#define UFILE_H
 
#include "unif01.h"


unif01_Gen * ufile_CreateReadText (char *fname, long nbuf);



void ufile_InitReadText (void);



unif01_Gen * ufile_CreateReadBin (char *fname, long nbuf);



void ufile_InitReadBin (void);


void ufile_DeleteReadText (unif01_Gen *);



void ufile_DeleteReadBin (unif01_Gen *);


void ufile_Gen2Bin (unif01_Gen *gen, char *fname, double n, int r, int s);

 
#endif
 

