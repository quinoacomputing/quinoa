/*************************************************************************\
 *
 * Package:        TestU01
 * File:           scatter.c
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

#include "util.h"
#include "chrono.h"
#include "mystr.h"

#include "scatter.h"
#include "unif01.h"

#include <string.h>
#include <stdio.h>



#define NUM_CHAR 100              /* Maximum length of file names */

#define LEN 250                   /* Maximum length of strings */



/*---------------------------- extern variables ---------------------------*/

long scatter_N;
int scatter_t;
lebool scatter_Over;
int scatter_x;
int scatter_y;
double scatter_L [scatter_MAXDIM + 1];
double scatter_H [scatter_MAXDIM + 1];
lebool scatter_Lacunary;
long scatter_LacI [scatter_MAXDIM + 1];
double scatter_Width = -1.0;
double scatter_Height = -1.0;
scatter_OutputType scatter_Output;
long scatter_Nkept;





/*---------------------------- module variables ---------------------------*/

static double V [scatter_MAXDIM + 1];  /* A point */

static char Nin [NUM_CHAR + 1] = {0};    /* Name of data file */

static char Nout1[NUM_CHAR + 1] = { 0 };
static char Nout2[NUM_CHAR + 1] = { 0 };
static char Nout3[NUM_CHAR + 1] = { 0 }; /* Names of results file */

static char S[NUM_CHAR + 1] = { 0 };
static char str[NUM_CHAR + 1] = { 0 }; /* Working strings */

static char Title[LEN + 1] = { 0 };

static int precision;             /* Number of decimals for the points */

static chrono_Chrono *chro;       /* Timer */






/*-------------------------------- functions ------------------------------*/

static unif01_Gen * scatter_ReadData (
   unif01_Gen *gen,
   char *F               /* Input file name without its extension .dat */
   )
/*
 * Reads the data in file <F>.dat for scatter_PlotUnif, according to the
 * format of the figure in the documentation scatter.tex.
 */
{
   FILE *fin;                     /* Input file */
   int i, j;
   double xa, xb;
   unif01_Gen *genL;

   strncpy (Nin, F, (size_t) NUM_CHAR - 5);
   /* Add the extension .dat to input data file */
   strcat (Nin, ".dat");
   fin = util_Fopen (Nin, "r");

   fgets (S, NUM_CHAR, fin);
   j = sscanf (S, " %ld", &scatter_N);
   util_Assert (j > 0, "scatter_ReadData:   on reading scatter_N");

   fgets (S, NUM_CHAR, fin);
   j = sscanf (S, " %d", &scatter_t);
   util_Assert (j > 0, "scatter_ReadData:   on reading scatter_t");
   util_Assert (scatter_t <= scatter_MAXDIM,
                "scatter_ReadData:   scatter_t > scatter_MAXDIM");
   util_Assert (scatter_t > 1, "scatter_ReadData:   scatter_t < 2");

   fgets (S, NUM_CHAR, fin);
   util_ReadBool (S, &scatter_Over);

   fgets (S, NUM_CHAR, fin);
   j = sscanf (S, " %d %d", &scatter_x, &scatter_y);
   util_Assert (j > 0,
                "scatter_ReadData:   on reading scatter_x or scatter_y");
   util_Assert (scatter_x <= scatter_t,
                "scatter_ReadData:  scatter_x > scatter_t");
   util_Assert (scatter_y <= scatter_t,
                "scatter_ReadData:  scatter_y > scatter_t");

   /* By default, bounds are [0, 1] */
   for (i = 1; i < scatter_t; i++) {
      scatter_L[i] = 0.0;
      scatter_H[i] = 1.0;
   }

   do {
      fgets (S, NUM_CHAR, fin);
      j = sscanf (S, " %d %lf %lf", &i, &xa, &xb);
      util_Assert (j > 0,
         "scatter_ReadData:   on reading r, scatter_L[r], scatter_H[r]");
      util_Assert (i <= scatter_t,
                   "scatter_ReadData:   r > scatter_t");
      scatter_L[i] = xa;
      scatter_H[i] = xb;
      util_Assert (scatter_L[i] >= 0.0,
         "scatter_ReadData:   scatter_L[r] < 0");
      util_Assert (scatter_H[i] <= 1.0,
                   "scatter_ReadData:   scatter_H[r] > 1");
      util_Assert (scatter_L[i] < scatter_H[i],
                   "scatter_ReadData:   scatter_H[r] <= scatter_L[r]");
   } while (i != scatter_t);

   fgets (S, NUM_CHAR, fin);
   j = sscanf (S, " %lf %lf", &scatter_Width, &scatter_Height);
   util_Assert (j > 0,
      "scatter_ReadData:   on reading scatter_Width, scatter_Height");

   fgets (S, NUM_CHAR, fin);
   j = sscanf (S, " %12s", str);
   if (!strcmp (str, "latex"))
      scatter_Output = scatter_latex;
   else if (!strcmp (str, "gnu_term"))
      scatter_Output = scatter_gnu_term;
   else if (!strcmp (str, "gnu_ps"))
      scatter_Output = scatter_gnu_ps;
   else {
      util_Error ("scatter_ReadData:   on reading scatter_Output");
   }

   fgets (S, NUM_CHAR, fin);
   j = sscanf (S, " %d", &precision);
   util_Assert (j > 0, "scatter_ReadData:   on reading Precision");

   fgets (S, NUM_CHAR, fin);
   util_ReadBool (S, &scatter_Lacunary);

   if (scatter_Lacunary) {
      for (i = 0; i < scatter_t; i++) {
         fgets (S, NUM_CHAR, fin);
         j = sscanf (S, " %ld", &scatter_LacI[i]);
         util_Assert (j > 0,
            "scatter_ReadData:   on reading scatter_LacI[]");
      }
      genL = unif01_CreateLacGen (gen, scatter_t, scatter_LacI);
   } else
      genL = gen;

   util_Fclose (fin);
   return genL;
}


/*=========================================================================*/

static unif01_Gen * scatter_ReadDataInterac (
   unif01_Gen *gen
   )
/*
 * Reads the data for scatter_PlotUnifInterac, interactively on the terminal
 */
{
   char rep;
   int i, j;
   lebool erreur;
   char format[12];
   unif01_Gen *genL;

   do {
      erreur = FALSE;
      printf ("What kind of output?\n"
              "latex:     (l)\n"
              "gnu_ps:    (p)\n"
              "gnu_term:  (t)\n");
      fgets (S, NUM_CHAR, stdin);
      j = sscanf (S, " %1s", str);

      rep = str[0];
      switch ((unsigned) rep) {
      case 'l':
      case 'L':
         scatter_Output = scatter_latex;
         break;
      case 'p':
      case 'P':
         scatter_Output = scatter_gnu_ps;
         break;
      case 't':
      case 'T':
         scatter_Output = scatter_gnu_term;
         break;
      default:
	 printf ("Please, answer with one letter amongst l, p, t.\n");
         erreur = TRUE;
         break;
      }
   } while (erreur);

   sprintf (S, "%1d", NUM_CHAR - 5);
   strcpy (format, " %");
   strcat (format, S);
   strcat (format, "s");
   printf ("Name of output file (without extension): ");
   fgets (S, NUM_CHAR, stdin);
   j = sscanf (S, format, Nin);
   util_Assert (j > 0, "scatter_ReadDataInterac");

   printf ("Number of points: ");
   fgets (S, NUM_CHAR, stdin);
   j = sscanf (S, " %ld", &scatter_N);
   util_Assert (j > 0, "scatter_ReadDataInterac");

   printf ("Number of dimensions: ");
   fgets (S, NUM_CHAR, stdin);
   j = sscanf (S, " %d", &scatter_t);
   util_Assert (j > 0, "scatter_ReadDataInterac");

   printf ("Overlapping:\n TRUE (t)\n FALSE (f)\n  ");
   fgets (S, NUM_CHAR, stdin);
   j = sscanf (S, " %1s", str);
   if (str[0] == 't')
      scatter_Over = TRUE;
   else
      scatter_Over = FALSE;

   printf ("Which dimension for the x-axis: ");
   fgets (S, NUM_CHAR, stdin);
   j = sscanf (S, " %d", &scatter_x);
   util_Assert (j > 0, "scatter_ReadDataInterac");

   printf ("Which dimension for the y-axis: ");
   fgets (S, NUM_CHAR, stdin);
   j = sscanf (S, " %d", &scatter_y);
   util_Assert (j > 0, "scatter_ReadDataInterac");

   for (i = 1; i <= scatter_t; i++) {
      printf ("Lower bound for x%1d: ", i);
      fgets (S, NUM_CHAR, stdin);
      j = sscanf (S, " %lf", &scatter_L[i]);
      util_Assert (j > 0, "scatter_ReadDataInterac");

      printf ("Upper bound for x%1d: ", i);
      fgets (S, NUM_CHAR, stdin);
      j = sscanf (S, " %lf", &scatter_H[i]);
      util_Assert (j > 0, "scatter_ReadDataInterac");
      util_Assert (scatter_L[i] >= 0.0,
                   "scatter_ReadDataInterac:   scatter_L[r] < 0");
      util_Assert (scatter_H[i] <= 1.0,
                   "scatter_ReadDataInterac:   scatter_H[r] > 1");
      util_Assert (scatter_L[i] < scatter_H[i],
                   "scatter_ReadDataInterac:   scatter_H[r] >= scatter_L[r]");
   }
   scatter_Height = 13.0;
   scatter_Width = 13.0;


   printf ("Lacunary:\n TRUE (t)\n FALSE (f)\n  ");
   fgets (S, NUM_CHAR, stdin);
   j = sscanf (S, " %1s", str);
   if (str[0] == 't')
      scatter_Lacunary = TRUE;
   else
      scatter_Lacunary = FALSE;

   if (scatter_Lacunary) {
      for (i = 0; i < scatter_t; i++) {
         printf ("Lacunary index %1d: ", i + 1);
         fgets (S, NUM_CHAR, stdin);
         j = sscanf (S, " %ld", &scatter_LacI[i]);
         util_Assert (j > 0,
             "scatter_ReadDataInterac:   on reading scatter_LacI[]");
      }
      genL = unif01_CreateLacGen (gen, scatter_t, scatter_LacI);
   } else
      genL = gen;

   printf ("Number of decimals of precision : ");
   fgets (S, NUM_CHAR, stdin);
   j = sscanf (S, " %d", &precision);
   util_Assert (j > 0, "scatter_ReadDataInterac:   on reading Precision");
   return genL;
}


/*=========================================================================*/

static lebool Retenu (void)
/*
 * Returns TRUE if the point with coordinates V[1],..., V[scatter_t] is
 * inside the bounds to be plotted; FALSE otherwise.
 */
{
   int j;
   for (j = 1; j <= scatter_t; j++) {
      if (V[j] < scatter_L[j] || V[j] > scatter_H[j])
         return FALSE;
   }
   return TRUE;
}


/*=========================================================================*/

static void HeadGraphTex (
   FILE *f               /* Latex output file */
   )
/*
 * Write the Latex commands before the points to be plotted.
 */
{
   fprintf (f, "\\documentclass [11pt]{article}\n"
               "\\begin {document}\n\n"
               "\\def\\fiverm {}%%\n"
               "\\input prepictex.tex \\input pictex.tex "
               "\\input postpictex.tex\n");
   fprintf (f, "\\begin{figure} \\centering \\beginpicture\n"
               "\\setcoordinatesystem units <%6.2fcm,%6.2fcm>\n",
               scatter_Width, scatter_Height);
   fprintf (f, "\\setplotarea x from 0 to 1, y from 0 to 1\n"
               "\\axis bottom\n"
               "  label $u_{n}$\n"
               "  ticks withvalues %8.4g %8.4g ",
               scatter_L[scatter_x],  scatter_H[scatter_x]);
   fprintf (f, " / at 0.0 1.0 / / \n"
               "\\axis left\n"
               "  label \\makebox[0pt]{$u_{n+%1d}$}\n", scatter_y - scatter_x);
   fprintf (f, "  ticks withvalues  %8.4g %8.4g ",
               scatter_L[scatter_y], scatter_H[scatter_y]);
   fprintf (f, " / at 0.0 1.0 / / \n"
               "\\axis top /  \\axis right /\n"
               "\\multiput {\\bf .} at\n");
}


/*-------------------------------------------------------------------------*/

static void BottomGraphTex (
   unif01_Gen *gen,
   FILE *f               /* Latex output file */
   )
/*
 * Write the Latex commands after the points have been written.
 */
{
   const double Epsilon = 1.0E-32;
   int j;
   char *p;
   size_t len;

   fprintf (f, "/ \\endpicture\n\n");
   /* fprintf (f, "\\parindent=0pt\n"); */
   fprintf (f,
           "\\def\\bornes#1#2#3 {$\\null\\ \\ \\ #2 < u_{n+#1} < #3$\\\\ }\n"
           "\\def\\bornez#1#2 {$\\null\\ \\ \\ #1 < u_{n} < #2$\\\\ }\n"
           "\\def\\stat#1#2 {\n"
           "Number of vectors generated: \\hbox to 1in {\\hfil #1.}\\\\\n"
           "Number of points plotted: \\hbox to 1in {\\hfil #2.}\\\\ }\n\n"
           "\\bigskip\\noindent {\\bf Generator:} \n");

   if ((p = strstr (gen->name, "Read"))) {
      /* This generator is from ufile. Write complete generator's name */
      strncpy (Title, gen->name, (size_t) LEN);
      p = strchr (gen->name, '\n');
      if (p)
         mystr_Subst (Title, "\n", "\n\n");
   } else {
      /* Remove initial seeds from gen name, but write filter if any */
      p = strchr (gen->name, ':');
      len = p - gen->name;
      strncpy (Title, gen->name, (size_t) len);
      Title[len] = '\0';
      p = strchr (gen->name, '\n');
      if (p) {
         strcat (Title, "\n");
         strcat (Title, p);
      }
   }

   /* Replace the _ in the generator name by \_ for Latex */
   mystr_Subst (Title, "_", "\\_");
   mystr_Subst (Title, "01_", "01\\_");
   fprintf (f, Title);

   fprintf (f, "\n\nHypercube in %1d dimensions.\\\\\n", scatter_t);
   fprintf (f, " Over = ");
   if (scatter_Over)
      fprintf (f, "TRUE");
   else
      fprintf (f, "FALSE");

   fprintf (f, "\\\\\n");
   for (j = 1; j <= scatter_t; j++) {
      if (scatter_L[j] > Epsilon || 1.0 - scatter_H[j] > Epsilon) {
	 if (j == 1) {
	    fprintf (f, "\\bornez {%9.4G}{%9.4G}\n",
                        scatter_L[j], scatter_H[j]);
	 } else {
	    fprintf (f, "\\bornes {%1d}{%9.4G}{%9.4G}\n",
                        j - 1, scatter_L[j], scatter_H[j]);
	 }
      }
   }
   /*
   if (scatter_Lacunary) {
      fprintf (f, "Lacunary indices = \\{");
      for (j = 0; j < scatter_t; j++) {
	if (j < scatter_t - 1)
           fprintf (f, " %ld,", scatter_LacI[j]);
        else
           fprintf (f, " %ld", scatter_LacI[j]);
      }
      fprintf (f, " \\}\\\\\n");
   }
   */
   fprintf (f, "\\stat {%10ld}{%10ld}\n", scatter_N, scatter_Nkept);
   fprintf (f, "Total CPU time : %12.2f", chrono_Val (chro, chrono_sec));
   fprintf (f, " seconds.\n"
               "\\end {figure}\n"
               "\\end {document}\n");
}


/*-------------------------------------------------------------------------*/

static void PutPointsTex (
   unif01_Gen *gen,
   FILE *f,              /* Latex output file */
   int Prec              /* Write the points with Prec decimals */
   )
{
   int j;
   long Npoints;                  /* Number of points */

   sprintf (S, "%%%1d", Prec + 5);
   sprintf (str, ".%1df", Prec);
   strcat (S, str);
   for (j = 1; j <= scatter_t; j++)
      V[j] = unif01_StripD (gen, 0);
   Npoints = 0;
   scatter_Nkept = 0;
   while (Npoints < scatter_N) {
      ++Npoints;
      if (Retenu ()) {
         ++scatter_Nkept;
         fprintf (f, S, (V[scatter_x] - scatter_L[scatter_x]) /
            (scatter_H[scatter_x] - scatter_L[scatter_x]));
         fprintf (f, S, (V[scatter_y] - scatter_L[scatter_y]) /
            (scatter_H[scatter_y] - scatter_L[scatter_y]));
         fprintf (f, "\n");
      }
      if (scatter_Over) {
         for (j = 1; j < scatter_t; j++)
            V[j] = V[j + 1];
         V[scatter_t] = unif01_StripD (gen, 0);
      } else {
         for (j = 1; j <= scatter_t; j++)
            V[j] = unif01_StripD (gen, 0);
      }
   }
}


/*=========================================================================*/

static void HeadGraphGnu (
   unif01_Gen *gen,
   char *F               /* File name without extension */
   )
/*
 * Creates the following output file names for gnuplot:
 * <F>.gnu contains the gnuplot commands to create the plot,
 * <F>.gnu.points contains the points to be plotted,
 * <F>.ps  contains the plot figure in PostScript format.
 * Write gnuplot commands for plot.
 */
{
   FILE *f;
   char *p, *q;
   size_t len;

   /* File of commands for gnuplot */
   strcpy (Nout1, F);
   strcat (Nout1, ".gnu");
   f = util_Fopen (Nout1, "w");

   /* File of points for gnuplot */
   strcpy (Nout2, Nout1);
   strcat (Nout2, ".points");

   fprintf (f, "set nokey\n" "set title \"");

   if ((p = strstr (gen->name, "Read"))) {
      /* This generator is from ufile. Write complete generator's name */
      strncpy (Title, gen->name, (size_t) LEN);
   } else {
      /* Remove initial seeds from gen name */
      p = strchr (gen->name, ':');
      len = p - gen->name;
      strncpy (Title, gen->name, (size_t) len);
      Title[len] = '\0';
   }
   /* Search for '\n' in title. If it is there, it will not be understood by 
      gnuplot. Process it specially to print the generator name */
   p = strchr (gen->name, '\n');
   if (p) {
      strncat (Title, p, (size_t) LEN);
      p = strchr (Title, '\n');
      q = Title;
      while (p) {
         *p = '\0';
         len = strlen (q);
         if (len > 0) {
            fprintf (f, q);
            fprintf (f, ";\\n");
         }
         p++;
         q = p;
         p = strchr (q, '\n');
      }
      fprintf (f, q);
   } else
      fprintf (f, Title);

   fprintf (f, ";\\n   N = %1ld", scatter_N);
   fprintf (f, "; t = %1d", scatter_t);
   if (scatter_Over)
      fprintf (f, "; Over");
   /*
   if (scatter_Lacunary) {
      int j;
      fprintf (f, "; LacI = {");
      for (j = 0; j < scatter_t; j++) {
         fprintf (f, " %ld", scatter_LacI[j]);
      }
      fprintf (f, " }");
   }
   */
   fprintf (f, "\"\nset xlabel \"u(n)\"\n"
               "set ylabel \"u(n+%1d)\"\n", scatter_y - scatter_x);
   fprintf (f, "set xrange [%4.2g:%4.2g]\n",
               scatter_L[scatter_x], scatter_H[scatter_x]);
   fprintf (f, "set yrange [%4.2g:%4.2g]\n",
               scatter_L[scatter_y], scatter_H[scatter_y]);
   fprintf (f, "set size square\n");
   if (scatter_Output == scatter_gnu_ps) {
      /* Set the size of the figure */
      /* fprintf (f, "set size %6.2f, %6.2f\n", Width/25.4, Height/17.8); */
      strcpy (Nout3, F);
      strcat (Nout3, ".ps");
      /* Postscript file for figure */
      fprintf (f, "set output \"");
      fprintf (f, Nout3);
      fprintf (f, "\"\nset term postscript");
   } else if (scatter_Output == scatter_gnu_term) {
      fprintf (f, "set output\n");
      fprintf (f, "set term x11");
   }
   fprintf (f, "\nplot \"");
   fprintf (f, Nout2);
   fprintf (f, "\"\n");
   if (scatter_Output == scatter_gnu_term) {
      fprintf (f, "pause -1  \"Hit return to continue \"\n");
   }
}


/*-------------------------------------------------------------------------*/

static void PutPointsGnu (
   unif01_Gen *gen,
   FILE *f,              /* Output file for the points in gnuplot format */
   int Prec              /* Write the points with Prec decimals */
   )
{
   int j;
   long Npoints;                  /* Number of points */

   sprintf (S, "%%%1d", Prec + 5);
   sprintf (str, ".%1df", Prec);
   strcat (S, str);
   for (j = 1; j <= scatter_t; j++)
      V[j] = unif01_StripD (gen, 0);
   Npoints = 0;
   scatter_Nkept = 0;
   while (Npoints < scatter_N) {
      ++Npoints;
      if (Retenu ()) {
         ++scatter_Nkept;
         fprintf (f, S, V[scatter_x]);
         fprintf (f, S, V[scatter_y]);
         fprintf (f, "\n");
      }
      if (scatter_Over) {
         for (j = 1; j < scatter_t; j++)
            V[j] = V[j + 1];
         V[scatter_t] = unif01_StripD (gen, 0);
      } else {
         for (j = 1; j <= scatter_t; j++)
            V[j] = unif01_StripD (gen, 0);
      }
   }
}


/*=========================================================================*/

static void Plot (unif01_Gen * gen, char *Nin, int Prec)
{
   FILE *f;

   if (scatter_Output == scatter_latex) {
      strcpy (Nout1, Nin);
      strcat (Nout1, ".tex");
      f = util_Fopen (Nout1, "w");
      HeadGraphTex (f);
      PutPointsTex (gen, f, Prec);
      BottomGraphTex (gen, f);
      util_Fclose (f);

   } else if (scatter_Output == scatter_gnu_ps ||
      scatter_Output == scatter_gnu_term) {
      HeadGraphGnu (gen, Nin);
      f = util_Fopen (Nout2, "w");
      PutPointsGnu (gen, f, Prec);
      util_Fclose (f);

   } else
      util_Error ("Plot:   scatter_Output has invalid value");
}


/*=========================================================================*/

void scatter_PlotUnif (unif01_Gen *gen, char *Nin)
{
   unif01_Gen *genL;
   genL = scatter_ReadData (gen, Nin);
   chro = chrono_Create ();
   Plot (genL, Nin, precision);
   chrono_Delete (chro);
}


/*=========================================================================*/

void scatter_PlotUnifInterac (unif01_Gen *gen)
{
   unif01_Gen *genL;
   genL = scatter_ReadDataInterac (gen);
   chro = chrono_Create ();
   Plot (genL, Nin, precision);
   chrono_Delete (chro);
}


/*=========================================================================*/

void scatter_PlotUnif1 (unif01_Gen * gen, long N, int Dim, lebool Over,
   int Proj[2], double Lower[], double Upper[], scatter_OutputType Output,
   int Prec, lebool Lac, long LacI[], char *Name)
{
   int j;
   unif01_Gen *genL;
   chro = chrono_Create ();
   scatter_N = N;
   scatter_t = Dim;
   scatter_Over = Over;
   scatter_x = Proj[0];
   scatter_y = Proj[1];
   for (j = 1; j <= scatter_t; j++) {
      scatter_L[j] = Lower[j - 1];
      scatter_H[j] = Upper[j - 1];
      util_Assert (scatter_L[j] >= 0.0, "scatter_PlotUnif1:   Lower[r] < 0");
      util_Assert (scatter_H[j] <= 1.0, "scatter_PlotUnif1:   Upper[r] > 1");
      util_Assert (scatter_L[j] < scatter_H[j],
                   "scatter_PlotUnif1:   Upper[r] <= Lower[r]");
   }
   if (scatter_Width <= 0.0)
      scatter_Width = 13.0;       /* cm */
   if (scatter_Height <= 0.0)
      scatter_Height = 13.0;      /* cm */
   scatter_Output = Output;
   scatter_Lacunary = Lac;
   if (scatter_Lacunary) {
      for (j = 0; j < scatter_t; j++)
         scatter_LacI[j] = LacI[j];
      genL = unif01_CreateLacGen (gen, scatter_t, scatter_LacI);
   } else
      genL = gen;
   strncpy (Nin, Name, (size_t) NUM_CHAR - 5);
   Plot (genL, Nin, Prec);
   chrono_Delete (chro);
}


/*=========================================================================*/
