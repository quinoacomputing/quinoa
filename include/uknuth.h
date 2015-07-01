
 
/* uknuth.h for ANSI C */

#ifndef UKNUTH_H
#define UKNUTH_H
 
#include "unif01.h"


unif01_Gen * uknuth_CreateRan_array1 (long s, long A[100]);



unif01_Gen * uknuth_CreateRan_array2 (long s, long A[100]);



unif01_Gen * uknuth_CreateRanf_array1 (long s, double B[100]);



unif01_Gen * uknuth_CreateRanf_array2 (long s, double B[100]);


void uknuth_DeleteRan_array1  (unif01_Gen *gen);
void uknuth_DeleteRan_array2  (unif01_Gen *gen);
void uknuth_DeleteRanf_array1 (unif01_Gen *gen);
void uknuth_DeleteRanf_array2 (unif01_Gen *gen);


 
#endif
 

