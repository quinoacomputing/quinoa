 
/* tables.h for ANSI C */
#ifndef TABLES_H
#define TABLES_H
 
#include "gdef.h"


typedef enum {
   tables_Plain,
   tables_Mathematica,
   tables_Matlab
   } tables_StyleType;


long ** tables_CreateMatrixL  (int M, int N);
unsigned long ** tables_CreateMatrixUL (int M, int N);
double ** tables_CreateMatrixD  (int M, int N);



void tables_DeleteMatrixL  (long *** T);
void tables_DeleteMatrixUL (unsigned long *** T);
void tables_DeleteMatrixD  (double *** T);



void tables_CopyTabL (long T1[], long T2[], int n1, int n2);
void tables_CopyTabD (double T1[], double T2[], int n1, int n2);



void tables_QuickSortL (long T[], int n1, int n2);
void tables_QuickSortD (double T[], int n1, int n2);

#ifdef USE_LONGLONG
   void tables_QuickSortLL (longlong T[], int n1, int n2);
   void tables_QuickSortULL (ulonglong T[], int n1, int n2);
#endif



void tables_WriteTabL (long V[], int n1, int n2, int k, int p, char Desc[]);

#ifdef USE_LONGLONG
   void tables_WriteTabLL (longlong V[], int n1, int n2, int k, int p,
                           char Desc[]);
   void tables_WriteTabULL (ulonglong V[], int n1, int n2, int k, int p,
                            char Desc[]);
#endif



void tables_WriteTabD (double V[], int n1, int n2, int k, int p1, int p2,
                       int p3, char Desc[]);



void tables_WriteMatrixD (double** Mat, int i1, int i2, int j1, int j2,
                          int w, int p, tables_StyleType style,
                          char Name[]);



void tables_WriteMatrixL (long** Mat, int i1, int i2, int j1, int j2, int w,
                          tables_StyleType style, char Name[]);



long tables_HashPrime (long n, double load);
 

#endif
 

