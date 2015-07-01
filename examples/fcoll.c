#include "fcong.h"
#include "ffam.h"
#include "fcho.h"
#include "fmultin.h"
#include "smultin.h"

int main (void)
{
   int NbDelta = 1;
   double ValDelta[] = { -1 };
   int t = 2;
   ffam_Fam *fam;
   smultin_Param *par;
   fmultin_Res *res;
   fcho_Cho *chon;
   fcho_Cho *chod;
   fcho_Cho2 *cho;

   fam = fcong_CreateLCG ("LCGGood.par", 10, 30, 1);
   par = smultin_CreateParam (NbDelta, ValDelta, smultin_GenerCellSerial, 2);
   res = fmultin_CreateRes (par);
   chon = fcho_CreateSampleSize (0.5, 1, 0, NULL, "n");
   chod = fmultin_CreatePer_DT (t, 1);
   cho = fcho_CreateCho2 (chon, chod);

   fmultin_Serial1 (fam, par, res, cho, 1, 0, t, TRUE, 21, 1, 5, 1);

   fcho_DeleteCho2 (cho);
   fmultin_DeletePer (chod);
   fcho_DeleteSampleSize (chon);
   fmultin_DeleteRes (res);
   smultin_DeleteParam (par);
   fcong_DeleteLCG (fam);
   return 0;
}
