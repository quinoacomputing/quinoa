#include "unif01.h"
#include "usoft.h"
#include "scatter.h"

int Proj[] = { 1, 3 };
long LacI[] = { 1, 2, 6};
double Lower[] = { 0.0, 0.0, 0.0 };
double Upper[] = { 0.0001, 0.5, 1.0 };

int main (void)
{
   unif01_Gen *gen;

   gen = usoft_CreateVisualBasic (12345);
   scatter_PlotUnif1 (gen, 10000000, 3, FALSE, Proj, Lower, Upper,
		      scatter_latex, 8, TRUE, LacI, "bone");
   usoft_DeleteGen (gen);
   return 0;
}
