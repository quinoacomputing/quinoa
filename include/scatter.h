
 
/* scatter.h  for ANSI C */
#ifndef SCATTER_H
#define SCATTER_H
 
#include "gdef.h"
#include "unif01.h"


typedef enum {
   scatter_latex,                 /* Latex format */
   scatter_gnu_ps,                /* gnuplot format for Postscript file */   
   scatter_gnu_term               /* Interactive gnuplot format */
   } scatter_OutputType;


#define scatter_MAXDIM 64


extern long scatter_N;             /* Number of points generated */
extern long scatter_Nkept;         /* Number of points kept for plot */
extern int scatter_t;              /* Dimension of points */
extern lebool scatter_Over;        /* = TRUE: overlapping points */
extern int scatter_x;
extern int scatter_y;              /* The 2 coordinates to plot */
extern double scatter_L [scatter_MAXDIM + 1];
extern double scatter_H [scatter_MAXDIM + 1];
                                   /* Lower and upper bounds for coordinates
                                      of points to plot */

extern lebool scatter_Lacunary;    /* = TRUE: lacunary case */
extern long scatter_LacI [scatter_MAXDIM + 1];  /* Lacunary indices */
extern double scatter_Width;
extern double scatter_Height;      /* Physical dimensions (in cm) of plot */
extern scatter_OutputType scatter_Output;       /* Kind of output */


void scatter_PlotUnif (unif01_Gen *gen, char *F);



void scatter_PlotUnif1 (unif01_Gen *gen, long N, int t, lebool Over,
   int Proj[2], double Lower[], double Upper[], scatter_OutputType Output,
   int Prec, lebool Lac, long LacI[], char *Name);



void scatter_PlotUnifInterac (unif01_Gen *gen);

 
#endif
 

