#include "unif01.h"
#include "ufile.h"
#include "scatter.h"

int main (void) 
{
   unif01_Gen *gen;
   gen = ufile_CreateReadText ("excel.pts", 100000);
   scatter_PlotUnif (gen, "excel");
   ufile_DeleteReadText (gen);
   return 0;
}
