 
/* uautomata.h for ANSI C */
#ifndef UAUTOMAT_H
#define UAUTOMAT_H
 
#include "unif01.h"

unif01_Gen * uautomata_CreateCA1 (int N, int S[ ], int r, int F[ ],
                                  int k, int ts, int cs, int rot);



unif01_Gen * uautomata_CreateCA90mp (int m, int S[]);


void uautomata_DeleteCA90mp (unif01_Gen *gen);



void uautomata_DeleteGen (unif01_Gen *gen);
 
#endif
 
