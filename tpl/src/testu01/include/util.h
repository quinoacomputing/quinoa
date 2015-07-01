 
/* util.h  for ANSI C */
#ifndef UTIL_H
#define UTIL_H
 
#include "gdef.h"
#include <stdio.h>
#include <stdlib.h>

 
#define util_Error(S) do { \
   printf ("\n\n******************************************\n"); \
   printf ("ERROR in file %s   on line  %d\n\n", __FILE__, __LINE__); \
   printf ("%s\n******************************************\n\n", S); \
   exit (EXIT_FAILURE); \
   } while (0)
 

 
#define util_Assert(Cond,S) if (!(Cond)) util_Error(S)
 

 
#define util_Warning(Cond,S) do { \
   if (Cond) { \
      printf ("*********  WARNING "); \
      printf ("in file  %s  on line  %d\n", __FILE__, __LINE__); \
      printf ("*********  %s\n", S); } \
   } while (0)
 

 
#define util_Max(x,y) (((x) > (y)) ? (x) : (y))
 

 
#define util_Min(x,y) (((x) < (y)) ? (x) : (y))
 


FILE * util_Fopen (const char *name, const char *mode);



int util_Fclose (FILE *stream);



int util_GetLine (FILE *file, char *Line, char c);



void util_ReadBool (char S[], lebool *x);



void util_WriteBool (lebool x, int d);



void * util_Malloc (size_t size);



void * util_Calloc (size_t dim, size_t size);



void * util_Realloc (void *ptr, size_t size);



void * util_Free (void *p);


 
#endif
 

