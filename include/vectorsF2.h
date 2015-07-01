 
#ifndef VECTORSF2_H
#define VECTORSF2_H
 
#include "gdef.h"


#define vectorsF2_WL      32


typedef struct{        
   int n;
   unsigned long *vect;
} BitVect;



void AllocBV ( BitVect *A,
               int l 
             );



void FreeBV ( BitVect *A 
            );



void BVCanonic ( BitVect *B, 
                 int t
               );



void AllOnes ( BitVect *B
             );



void Mask ( BitVect *B,
            int t
          );
void InvMask ( BitVect *B,
               int t
             );



lebool BVisZero ( BitVect *A
                );



int ValBitBV ( BitVect *A, 
               int noBit
             );                     



void PutBitBV ( BitVect *A,
                int noBit, 
                int valBit
              );           



void PutBVToZero ( BitVect *A
                 );



void CopyBV ( BitVect *A, 
              BitVect *B
            );
   


void CopyBVPart ( BitVect *A,
                  BitVect *B,
                  int l
                );



lebool CompareBV ( BitVect *A, 
                    BitVect *B
                  ); 



void XORBV ( BitVect *A, 
             BitVect *B, 
             BitVect *C
           );



void XOR2BV ( BitVect *A, 
              BitVect *B, 
              BitVect *C, 
              BitVect *D
            );



void ANDBV ( BitVect *A, 
             BitVect *B, 
             BitVect *C
           );



void ANDBVSelf ( BitVect *A, 
                 BitVect *B
               );



void ANDBVMask ( BitVect *A,
                 BitVect *B,
                 int t
               );



void ANDBVInvMask ( BitVect *A,
                    BitVect *B,
                    int t
                  );



void XORBVSelf ( BitVect *A, 
                 BitVect *B
               );



void BVLShift ( BitVect *R, 
                BitVect *A, 
                int n 
              );    



void BVRShift ( BitVect *R, 
                BitVect *A, 
                int n 
              );



void BVLShiftSelf ( BitVect *R, 
                    int n 
                  );



void BVLS1Self ( BitVect *R 
               );



void BVRShiftSelf ( BitVect *R, 
                    int n 
                  );



void InverseBV ( BitVect *A
               );



void DispBitVect ( BitVect *A, 
                   int l,
                   int mathematica
                 );



void RandVect ( BitVect *v
              );



lebool VerifBitsCommuns ( BitVect *ds1, 
                           BitVect *ds2 
                         );



typedef struct{
  BitVect **lignes;
  int nblignes;
  int t;
  int l;
} Matrix;



void AllocMat ( Matrix* m,
                int nblines,
                int l,
                int t
              );



void FreeMat ( Matrix *m
             );



void CopyMat ( Matrix *m, 
               Matrix *ms, 
               int nl, 
               int t 
             );



void CopyNTupleMat ( Matrix *m, 
                     Matrix *ms, 
                     int nl,
                     int *colonnes, 
                     int t 
                   );



lebool Diag ( Matrix *m, 
               int kg,                                
               int t, 
               int l, 
               int *gr 
             );



int GaussianElimination ( Matrix *m, 
                          int nblignes,
                          int l,
                          int t
                        );



int CompleteElimination ( Matrix *m, 
                          int nblignes,
                          int l,
                          int t
                        );



int SpecialGaussianElimination ( Matrix *m, 
                                 int nblignes, 
                                 int l, 
                                 int t, 
                                 int *indices
                               );



void MultMatrixByBV ( BitVect *A,
                     Matrix *m, 
                     BitVect *B
                   );



void TransposeMatrices ( Matrix *T, 
                         Matrix *M, 
                         int nblines, 
                         int t, 
                         int l
                       );



void ExchangeVect ( Matrix *m, 
                    int i, 
                    int j 
                  );



void XorVect ( Matrix *m,    
               int r, 
               int s, 
               int min, 
               int max
             );



void DispMat ( Matrix *m, 
               int t, 
               int l,
               int nblines,
               lebool mathematica
            );
 


lebool InverseMatrix ( Matrix *MInv, 
                        Matrix *M
                      );



void MultMatrixByMatrix ( Matrix *A, 
                          Matrix *B,
                          Matrix *C
                        );



void MatrixTwoPow ( Matrix *A, 
                    Matrix *B,
                    unsigned int e
                  );


#ifdef USE_LONGLONG
   void MatrixPow (Matrix * A, Matrix * B, longlong e);
#else
   void MatrixPow (Matrix * A, Matrix * B, long e);
#endif

 
#endif
 

