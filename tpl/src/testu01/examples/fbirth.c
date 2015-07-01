#include "fcong.h"
#include "ffam.h"
#include "fcho.h"
#include "fmarsa.h"

int main (void)
{
   ffam_Fam *fam;
   fcho_Cho *chon;
   fcho_Cho *chod;
   fcho_Cho2 *cho;

   fam = fcong_CreateLCGPow2 (NULL, 10, 30, 1);
   chon = fcho_CreateSampleSize (1.0/3.0, 1, 0, NULL, "n");
   chod = fmarsa_CreateBirthEC (1, 2, 1.0);
   cho = fcho_CreateCho2 (chon, chod);
   fmarsa_BirthdayS1 (fam, NULL, cho, 1, 0, 2, 1, 21, 1, 5, 1);
   fcho_DeleteCho2 (cho);
   fmarsa_DeleteBirthEC (chod);
   fcho_DeleteSampleSize (chon);
   fcong_DeleteLCGPow2 (fam);
   return 0;
}
