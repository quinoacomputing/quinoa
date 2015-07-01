 
/* fspectral.h  for ANSI C */
#ifndef FSPECTRAL_H
#define FSPECTRAL_H
 
#include "ffam.h"
#include "fres.h"
#include "fcho.h"


extern long fspectral_MaxN;

void fspectral_Fourier3 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         int k, int r, int s,
                         int Nr, int j1, int j2, int jstep);

 
#endif
 

