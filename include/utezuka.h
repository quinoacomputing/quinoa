
 
/* utezuka.h for ANSI C */

#ifndef UTEZUKA_H
#define UTEZUKA_H
 
#include "unif01.h"


unif01_Gen * utezuka_CreateTezLec91 (unsigned int Y1, unsigned int Y2);



unif01_Gen * utezuka_CreateTez95 (unsigned int Y1, unsigned int Y2,
                                  unsigned int Y3);



unif01_Gen * utezuka_CreateTezMRG95 (unsigned int Y1[5],
                                     unsigned int Y2[7]);


void utezuka_DeleteGen (unif01_Gen *gen);


 
#endif
 

