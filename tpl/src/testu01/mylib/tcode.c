/*************************************************************************\
 *
 * Package:        MyLib
 * File:           tcode.c
 * Environment:    ANSI C
 *
 * Copyright (c) 2002 Pierre L'Ecuyer, DIRO, Université de Montréal.
 * e-mail: lecuyer@iro.umontreal.ca
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted without a fee for private, research,
 * academic, or other non-commercial purposes.
 * Any use of this software in a commercial environment requires a
 * written licence from the copyright owner.
 *
 * Any changes made to this package must be clearly identified as such.
 *
 * In scientific publications which used this software, a reference to it
 * would be appreciated.
 *
 * Redistributions of source code must retain this copyright notice
 * and the following disclaimer.
 *
 * THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
\*************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>





/*========================== constants =================================*/

#define MaxChar 255                  /* max. length of a line */

#define Mess1 "\nUsage:   tcode  <FileIn>  <FileOut>\n\n"
#define Mess2 "\nERROR:   The 2 files must be different\n\n"




/*============================= variables ==============================*/

static char Line[MaxChar + 1] = {0}, /* One line of input */
   FIn[MaxChar + 1] = {0},           /* Input file name */
   FOut[MaxChar + 1] = {0};          /* Output file name */

static FILE *fin,                    /* Input file */
  *fout;                             /* Output file */

static size_t L1, L2, L3, L4, L5, L6, L7, L8; /* Lengths of strings */




/**************************************************************************/

static void Init (int argc, char *argv[])
{
   if (argc != 3) {
      printf (Mess1);
      exit (EXIT_FAILURE);
   }
   sprintf (FIn, "%.*s", MaxChar, argv[1]);
   sprintf (FOut, "%.*s", MaxChar, argv[2]);
   if (!strcmp (FIn, FOut)) {
      printf (Mess2);
      exit (EXIT_FAILURE);
   }

   errno = 0;
   fin = fopen (FIn, "r");
   if (fin == NULL) {
      printf ("\nOpening of %s failed: %s\n\n", FIn, strerror (errno));
      exit (EXIT_FAILURE);
   }
   errno = 0;
   fout = fopen (FOut, "w");
   if (fout == NULL) {
      printf ("\nOpening of %s failed: %s\n\n", FOut, strerror (errno));
      exit (EXIT_FAILURE);
   }

   L1 = strlen ("\\code");
   L2 = strlen ("\\endcode");
   L3 = strlen ("\\hide");
   L4 = strlen ("\\endhide");
   L5 = strlen ("\\iffalse");
   L6 = strlen ("\\fi");
   L7 = strlen ("\\smallc");
   L8 = strlen ("\\smallcode");
}

/**************************************************************************/

static void SearchRep (char *A, const char *sub, unsigned int L)
{
   /* Search for TEX command sub, of length L, in line A */
   char *p, *q;

   p = strstr (A, sub);              /* search for sub in A */
   if (p) {                          /* if found */
      if (isalpha ((int) *(p + L)))  /* is it only part of a command? */
         SearchRep (A + L, sub, L);  /* yes? search rest of line */
      else {
         q = p + L;                  /* no? delete */
         *p++ = ' ';
         while ((*p++ = *q++))
            ;
      }
   }
}

/**************************************************************************/

static void ProcessCode (char *A)
{

   SearchRep (A, "\\hide", L3);
   SearchRep (A, "\\endhide", L4);
   SearchRep (A, "\\iffalse", L5);
   SearchRep (A, "\\fi", L6);
   SearchRep (A, "\\smallcode", L8);
   SearchRep (A, "\\smallc", L7);

   fprintf (fout, "%s", A);
}

/**************************************************************************/

static int ProcessLine (char *line)
{
   /* Process a line (or part of a line) of valid code. If at the end of the
      line, we are still in a region of valid code, then return 1; otherwise,
      return 0. */

   char *p, *rest;
   int code;

   p = strstr (line, "\\endcode");   /* search for "\endcode" */
   if (p) {                          /* if found */
      *p = '\0';                     /* code ends here */
      rest = p + L2;                 /* step over "\endcode" */
      p = strstr (rest, "\\code");   /* search for "\code" */
      if (p) {
         p += L1;                    /* step over "\code" */
         rest = calloc (MaxChar, sizeof (char));
         strncpy (rest, p, MaxChar);
         ProcessCode (line);
         code = ProcessLine (rest);
         free (rest);
         return code;
      } else {
         ProcessCode (line);
         return 0;
      }
   } else {
      ProcessCode (line);
      return 1;
   }
}


/**************************************************************************/

int main (int argc, char *argv[])
{
   char *p, *q;
   int isCode = 0;         /* If isCode == TRUE, we are in a region of
                              valid code; otherwise not. */

   Init (argc, argv);

   while (NULL != fgets (Line, MaxChar, fin)) { /* Not EOF and no error */
      if (isCode) {
         isCode = ProcessLine (Line);
      } else {
         /* search for "\def\code" and drop that line: it is not valid code
            but the definition of the TEX command \code */
         if (strstr (Line, "\\def\\code"))
            ;                             
         else if ((p = strstr (Line, "\\code"))) {
            /* search for "\code". If "\code" is found on a line with a %
               before it, then it is a TEX comment and we do not consider
               it as starting a region of valid code; */
            *p = '\0';
            q = strchr (Line, '%');
            if (NULL == q)
               /* otherwise, it is valid code: process rest of line */
               isCode = ProcessLine (p + L1);
         }
      }
   }

   fprintf (fout, "\n");
   fclose (fout);
   fclose (fin);
   return 0;
}
