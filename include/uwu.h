 
/* uwu.h for ANSI C */

#ifndef UWU_H
#define UWU_H
 
#include "gdef.h"
#include "unif01.h"


#ifdef USE_LONGLONG
   unif01_Gen * uwu_CreateLCGWu61a (longlong s);



   unif01_Gen * uwu_CreateLCGWu61b (longlong s);

#endif



void uwu_DeleteGen (unif01_Gen *gen);


 
#endif
 

