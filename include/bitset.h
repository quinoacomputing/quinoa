 
/* bitset.h  for ANSI C */
#ifndef BITSET_H
#define BITSET_H
#include "gdef.h"
 

extern unsigned long bitset_maskUL[];



extern unsigned long bitset_MASK[];


typedef unsigned long bitset_BitSet;

 
#define bitset_SetBit(S, b) ((S) |= (bitset_maskUL[b]))
 

 
#define bitset_ClearBit(S, b) ((S) &= ~(bitset_maskUL[b]))
 

 
#define bitset_FlipBit(S, b) ((S) ^= (bitset_maskUL[b]))
 

 
#define bitset_TestBit(S, b)  ((S) & (bitset_maskUL[b]) ? 1 : 0)
 

 
#define bitset_RotateLeft(S, t, r) do { \
   unsigned long v853 = (S) >> ((t) - (r)); \
   (S) <<= (r);   (S) |= v853;   (S) &= bitset_MASK[t]; \
   } while (0)
 

 
#define bitset_RotateRight(S, t, r) do { \
   unsigned long v972 = (S) >> (r); \
   (S) <<= ((t) - (r));   (S) |= v972;   (S) &= bitset_MASK[t]; \
   } while (0)
 


bitset_BitSet bitset_Reverse (bitset_BitSet Z, int s);



void bitset_WriteSet (char *desc, bitset_BitSet Z, int s);

 
#endif
 

