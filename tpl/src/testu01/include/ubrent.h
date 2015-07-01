 
/* ubrent.h for ANSI C */
#ifndef UBRENT_H
#define UBRENT_H
 
#include "gdef.h"
#include "unif01.h"


unif01_Gen* ubrent_CreateXorgen32 (int r, int s, int a, int b, int c, int d,
                                   lebool hasWeyl, unsigned int seed);


unif01_Gen* ubrent_CreateXor4096s (unsigned int seed);



#ifdef USE_LONGLONG

unif01_Gen* ubrent_CreateXorgen64 (int r, int s, int a, int b, int c, int d,
                                   lebool hasWeyl, ulonglong seed);



unif01_Gen* ubrent_CreateXor4096l (ulonglong seed);



unif01_Gen* ubrent_CreateXor4096d (ulonglong seed);


#endif


unif01_Gen* ubrent_CreateXor4096i (unsigned long seed);



unif01_Gen* ubrent_CreateXor4096r (unsigned long seed);


void ubrent_DeleteXorgen32 (unif01_Gen *);
void ubrent_DeleteXor4096s (unif01_Gen *);
void ubrent_DeleteXor4096i (unif01_Gen *);
void ubrent_DeleteXor4096r (unif01_Gen *);

#ifdef USE_LONGLONG
   void ubrent_DeleteXorgen64 (unif01_Gen *);
   void ubrent_DeleteXor4096l (unif01_Gen *);
   void ubrent_DeleteXor4096d (unif01_Gen *);
#endif
 
#endif
 
