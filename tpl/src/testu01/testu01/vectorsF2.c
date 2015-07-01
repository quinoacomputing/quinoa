#include "vectorsF2.h"
#include <stdio.h>
#include <stdlib.h>

#define WL vectorsF2_WL

#define MC 0x80000000UL    /* permet de diagonaliser la matrice dans Diag() */

unsigned long MMC[WL] =
   { MC, MC >> 1, MC >> 2, MC >> 3, MC >> 4, MC >> 5, MC >> 6, MC >> 7,
   MC >> 8, MC >> 9, MC >> 10,
   MC >> 11, MC >> 12, MC >> 13, MC >> 14, MC >> 15, MC >> 16, MC >> 17,
   MC >> 18, MC >> 19, MC >> 20,
   MC >> 21, MC >> 22, MC >> 23, MC >> 24, MC >> 25, MC >> 26, MC >> 27,
   MC >> 28, MC >> 29, MC >> 30, MC >> 31
};

lebool InverseMatrix (Matrix * InvM, Matrix * M)
{

   Matrix Temp;
   int j, rang;
   if (M->nblignes != M->l) {
      printf ("Matrix M is not square!\n");
      exit (1);
   }
   AllocMat (&Temp, M->nblignes, M->l, 2);
   for (j = 0; j < M->l; j++)
      CopyBV (&(Temp.lignes[j][0]), &(M->lignes[j][0]));
   for (j = 0; j < M->l; j++) {
      BVCanonic (&(Temp.lignes[j][1]), j);
   }
   /* DispMat(&Temp,2,M->l,M->nblignes,0); */
   rang = CompleteElimination (&Temp, M->nblignes, M->l, 2);
   /* DispMat(&Temp,2,M->l,M->nblignes,0); */
   /* printf("rang=%d",rang); */
   for (j = 0; j < M->l; j++)
      CopyBV (&(InvM->lignes[j][0]), &(Temp.lignes[j][1]));
   return (rang == M->l);
   FreeMat (&Temp);
}


/* ********************************************************************** */
/* lebool Diag( Matrix m, int kg,				          */
/*               int t, int l, int *gr )                                  */
/* Evalue si la matrice de travail m sur kg lignes est de plein rang en   */
/* la diagonalisant.  On procede sur t BitVect en considerant les l       */
/* premiers bits de chacun.  La fonction retourne TRUE si la matrice m    */
/* est de plein rang t*l et *gr est inchange. La fonction retourne        */
/* FALSE sinon et *gr prend pour valeur le numero du BitVect ou il y a    */
/* eu echec moins un ( = dimension pour laquelle on a resolution l ).     */
/* ********************************************************************** */
lebool Diag (Matrix * m, int kg, int t, int l, int *gr)
{
   int i, j, cl, rang;

   rang = 0;

   /* On diagonalise la matrice des entrees (i,j) sur l bits */
   /* avec 0 <= i < kg et 0 <= j < t .  */

   for (j = 0; j < t; j++) {
      cl = 1;
      while (cl <= l) {
         /* On cherche dans la j-eme colonne, commencant a la */
         /* rang-eme ligne, la premiere entree dont le bit le plus */
         /* significatif est non nul.  Bref, on cherche un pivot.  */
         i = rang;

         while ((i < kg)
            && (m->lignes[i][j].vect[(cl - 1) / WL] < MMC[(cl - 1) % WL]))
            i++;
         if (i < kg) {            /* pivot trouve ... */
            ExchangeVect (m, rang, i);
            for (i = rang + 1; i < kg; i++) {
               if (m->lignes[i][j].vect[(cl - 1) / WL] & MMC[(cl - 1) % WL])
                  XorVect (m, i, rang, j, m->t);
            }
            rang++;
         } else {                 /* pas de pivot trouve ... ==> pas de plein
                                     rang ... */

            *gr = j;              /* no de groupe ou il y a echec moins un */
            return FALSE;         /* c'est j car on indexe a partir de 0 .  */
         }
         cl++;
      }
   }
   return TRUE;                   /* on a trouve tous les pivots ==> plein
                                     rang ! */
}


int CompleteElimination (Matrix * m, int nblignes, int l, int t)
{
   int i, j, cl, rang;

   rang = 0;

   j = 0;
   while (j < t) {
      cl = 0;
      while (cl < l) {
         /* On cherche dans la j-eme colonne, commencant a la */
         /* rang-eme ligne, la premiere entree dont le bit le plus */
         /* significatif est non nul.  Bref, on cherche un pivot.  */
         i = rang;
         while ((i < nblignes)
            && ((((m->lignes)[i])[j]).vect[(cl) / WL] < MMC[(cl) % WL]))
            i++;
         if (i < nblignes) {      /* pivot trouve ... */
            ExchangeVect (m, rang, i);
            for (i = 0; i < nblignes; i++)
               if (i != rang)
                  if ((((m->lignes)[i])[j]).vect[(cl) / WL] & MMC[(cl) % WL])
                     XorVect (m, i, rang, 0, m->t);

            rang++;
            if (rang == nblignes)
               return rang;
         } else
            return rang;
         cl++;
      }
      j++;
   }
   return rang;
}


int GaussianElimination (Matrix * m, int nblignes, int l, int t)
{
   int i, j, cl, rang;

   rang = 0;

   j = 0;
   while (j < t) {
      cl = 0;
      while (cl < l) {
         /* On cherche dans la j-eme colonne, commencant a la */
         /* rang-eme ligne, la premiere entree dont le bit le plus */
         /* significatif est non nul.  Bref, on cherche un pivot.  */
         i = rang;
         while ((i < nblignes)
            && ((((m->lignes)[i])[j]).vect[(cl) / WL] < MMC[(cl) % WL]))
            i++;
         if (i < nblignes) {      /* pivot trouve ... */
            ExchangeVect (m, rang, i);
            for (i = rang + 1; i < nblignes; i++)
               if ((((m->lignes)[i])[j]).vect[(cl) / WL] & MMC[(cl) % WL])
                  XorVect (m, i, rang, 0, m->t);

            rang++;
            if (rang == nblignes)
               return rang;
         }
         cl++;
      }
      j++;
   }
   return rang;
}


int SpecialGaussianElimination (Matrix * m, int nblignes, int l, int t,
   int *indices)
{
   int i, j, cl, rang;

   rang = 0;

   j = 0;
   while (j < t) {
      cl = 0;
      while (cl < l) {
         /* On cherche dans la j-eme colonne, commencant a la */
         /* rang-eme ligne, la premiere entree dont le bit le plus */
         /* significatif est non nul.  Bref, on cherche un pivot.  */
         i = rang;
         while ((i < nblignes)
            && ((((m->lignes)[i])[indices[j]]).vect[(cl) / WL] <
               MMC[(cl) % WL]))
            i++;
         if (i < nblignes) {      /* pivot trouve ... */
            ExchangeVect (m, rang, i);
            for (i = rang + 1; i < nblignes; i++)
               if ((((m->lignes)[i])[indices[j]]).vect[(cl) /
                     WL] & MMC[(cl) % WL])
                  XorVect (m, i, rang, 0, m->t);

            rang++;
            if (rang == nblignes)
               return rang;
         }
         cl++;
      }
      j++;
   }
   return rang;
}


void MultMatrixByBV (BitVect * A, Matrix * M, BitVect * B)
{
   int i, j, res;
   if (M->l < B->n * WL) {
      printf ("Error in MultMatrixByBV(): sizes do not match\n");
      exit (1);
   }
   if (A->n * WL < M->nblignes) {
      printf ("Error in MultMatrixByBV(): sizes do not match\n");
      exit (1);
   }
   if (M->t != 1) {
      printf ("Error in MultMatrixByBV(): Not implemented for M->t > 1\n");
      exit (1);
   }
   PutBVToZero (A);
   for (i = 0; i < M->nblignes; i++) {
      res = 0;
      for (j = 0; j < M->l; j++)
         res += ValBitBV (&(M->lignes[i][0]), j) * ValBitBV (B, j);
      res %= 2;
      PutBitBV (A, i, res);
   }
}


/* ********************************************************************** */
/* int GetBitBV (BitVect A, int noBit)			  */
/* Fonction qui permet de prendre la valeur du noBit-ieme bit             */
/* (l'indexation commence au bit 0)					  */
/* ********************************************************************** */
int ValBitBV (BitVect * A, int noBit)
{
   int k;
   unsigned long mask;
   k = noBit / WL;
   mask = 0x80000000UL >> (noBit - k * WL);
   if (A->vect[k] & mask)
      return 1;
   else
      return 0;
}


/* ********************************************************************** */
/* void SetBitBV (BitVect A, int noBit, int valBit)		  */
/* Fonction qui permet de rendre la valeur du noBit-ieme bit `a valBit     */
/* (l'indexation commence au bit 0) (valBit=1 ou valBit=0)		  */
/* ********************************************************************** */
void PutBitBV (BitVect * A, int noBit, int valBit)
{
   int k;
   unsigned long mask;

   k = noBit / WL;
   if (valBit == 1) {
      mask = 0x80000000UL >> (noBit - k * WL);
      A->vect[k] |= mask;
   } else {
      mask = 0xffffffffUL ^ (0x80000000UL >> (noBit - k * WL));
      A->vect[k] &= mask;
   }
}


/* ********************************************************************** */
/* SetBVToZero( BitVect A)        					  */
/* Initialise le vecteur de bit a zero       				  */
/* ********************************************************************** */
void PutBVToZero (BitVect * A)
{
   int i;
   for (i = 0; i < A->n; i++)
      A->vect[i] = 0UL;
}


/* ********************************************************************** */
/* void CopyBV(BitVect A, BitVect B)       				  */
/* Copie le contenu de B dans A (A=B)       				  */
/* ********************************************************************** */
void CopyBV (BitVect * A, BitVect * B)
{
   int i;

   if (A->n != B->n) {
      printf
    ("Error in CopyBV(): vectors of different dimensions! (%d and %d bits)\n",
         A->n * WL, B->n * WL);
      exit (1);
   }

   if (B->n == 0) {
      printf ("Nothing to copy!\n");
      exit (1);
   }
   for (i = 0; i < B->n; i++)
      A->vect[i] = B->vect[i];
}
void CopyBVPart (BitVect * A, BitVect * B, int l)
{

   int i, n;

   n = (l - 1) / WL + 1;

   if (A->n < n) {
      printf ("Error in CopyBVPart() : The vector A is not large enough!\n");
      exit (1);
   }
   if (B->n == 0) {
      printf ("Nothing to copy!\n");
      exit (1);
   }

   for (i = 0; i < n; i++)
      A->vect[i] = B->vect[i];

   if (l % WL) {
      BitVect m;
      AllocBV (&m, A->n * WL);
      Mask (&m, l);
      ANDBVSelf (A, &m);
      FreeBV (&m);
   }
}


/* ********************************************************************** */
/* void EgalBV(BitVect A, BitVect B)       				  */
/* Compare le contenu de B avec celui de A.  Retourne TRUE si les deux    */
/* contiennent la meme information.       				  */
/* ********************************************************************** */
lebool CompareBV (BitVect * A, BitVect * B)
{

   int i;

   if (A->n != B->n) {
      printf ("Error in EgalBV(): Vectors of different sizes\n");
      exit (1);
   }

   for (i = 0; i < A->n; i++)
      if (A->vect[i] != B->vect[i])
         return FALSE;
   return TRUE;
}


lebool BVisZero (BitVect * A)
{
   int j = 0;
   while (j < A->n)
      if (A->vect[j++] != 0UL)
         return FALSE;
   return TRUE;
}


/* ********************************************************************** */
/* void XORBV(BitVect A, BitVect B, BitVect C)      			  */
/* Cette fonction effectue A = B ^ C       				  */
/* ********************************************************************** */
void XORBV (BitVect * A, BitVect * B, BitVect * C)
{

   int i;

   if ((A->n != B->n) || (B->n != C->n)) {
      printf ("Error in XORBV(): Vectors of different sizes\n");
      exit (1);
   }

   for (i = 0; i < B->n; i++)
      A->vect[i] = B->vect[i] ^ C->vect[i];
}


/* ********************************************************************** */
/* void XOR2BV(BitVect A, BitVect B, BitVect C, BitVect D)    		  */
/* Cette fonction effectue A = B ^ C ^ D      				  */
/* ********************************************************************** */
void XOR2BV (BitVect * A, BitVect * B, BitVect * C, BitVect * D)
{

   int i;

   if ((A->n != B->n) || (B->n != C->n) || (C->n != D->n)) {
      printf ("Error in XOR2BV(): Vectors of different sizes\n");
      exit (1);
   }

   for (i = 0; i < B->n; i++)
      A->vect[i] = B->vect[i] ^ C->vect[i] ^ D->vect[i];
}


/* ********************************************************************** */
/* void ANDBV(BitVect A, BitVect B, BitVect C)      			  */
/* Cette fonction effectue A = B & C       				  */
/* ********************************************************************** */
void ANDBV (BitVect * A, BitVect * B, BitVect * C)
{

   int i;

   if ((A->n != B->n) || (B->n != C->n)) {
      printf ("Error in ANDBV(): Vectors of different sizes\n");
      exit (1);
   }

   for (i = 0; i < B->n; i++)
      A->vect[i] = B->vect[i] & C->vect[i];
}


void ANDBVMask (BitVect * A, BitVect * B, int t)
{
   int n, m, j;
   if (A->n != B->n) {
      printf ("Error in ANDBVMask(): Vectors of different sizes\n");
      exit (1);
   }

   if (t > B->n * WL)
      CopyBV (A, B);
   else if (t == 0)
      PutBVToZero (A);
   else {
      n = t / WL;
      m = t - n * WL;
      for (j = 0; j < n; j++) {
         A->vect[j] = B->vect[j];

      }
      if (m != 0) {
         A->vect[j] = B->vect[j] & (0xffffffffUL << (WL - m));
         j++;
      }
      /* printf("n=%d j=%d %d m=%d ",n,j,A->n,m); */
      for (; j < A->n; j++) {
         A->vect[j] = 0UL;
      }
   }
}


void ANDBVInvMask (BitVect * A, BitVect * B, int t)
{
   int n, m, j;
   if (A->n != B->n) {
      printf ("Error in ANDBV(): Vectors of different sizes\n");
      exit (1);
   }
   if (t > B->n * WL)
      PutBVToZero (A);
   else if (t == 0)
      CopyBV (A, B);
   else {
      n = t / WL;
      m = t - n * WL;
      for (j = 0; j < n; j++)
         A->vect[j] = 0UL;
      if (m == 0)
         A->vect[j] = B->vect[j];
      else {
         A->vect[j] = B->vect[j] & (0xffffffffUL >> m);

      }
      j++;
      for (; j < A->n; j++)
         A->vect[j] = B->vect[j];
   }
}




/* ********************************************************************** */
/* void ANDBVSelf(BitVect A, BitVect B)       				  */
/* Cette fonction effectue A &=B     					  */
/* ********************************************************************** */
void ANDBVSelf (BitVect * A, BitVect * B)
{
   int i;

   if ((A->n != B->n)) {
      printf ("Error in ANDBVSelf(): Vectors of different sizes\n");
      exit (1);
   }

   for (i = 0; i < B->n; i++)
      A->vect[i] &= B->vect[i];
}


/* ********************************************************************** */
/* void XORBVSelf(BitVect A, BitVect B)     				  */
/* Cette fonction effectue A ^=B    					  */
/* ********************************************************************** */
void XORBVSelf (BitVect * A, BitVect * B)
{
   int i;

   if ((A->n != B->n)) {
      printf ("Error in XORBVSelf(): Vectors of different sizes\n");
      exit (1);
   }

   for (i = 0; i < B->n; i++)
      A->vect[i] ^= B->vect[i];
}


/* ********************************************************************** */
/* void BVLShift( BitVect R, BitVect A, int n )                      */
/* Effectue : R = A << n ;                                                */
/* ********************************************************************** */
void BVLShift (BitVect * R, BitVect * A, int n)
{
   int i;
   int WLmn;
   unsigned long temp;

   if ((R->n != A->n)) {
      printf ("Error in BVLShift(): Vectors of different sizes\n");
      exit (1);
   }

   for (i = 0; i < A->n; i++)
      R->vect[i] = A->vect[i];
   while (n >= 32) {
      for (i = 1; i < A->n; i++)
         R->vect[i - 1] = R->vect[i];
      R->vect[A->n - 1] = 0UL;
      n -= 32;
   }
   if (n > 0) {
      WLmn = WL - n;
      R->vect[0] <<= n;
      for (i = 1; i < A->n; i++) {
         temp = R->vect[i] >> WLmn;
         R->vect[i - 1] |= temp;
         R->vect[i] <<= n;
      }
   }

}


/* ********************************************************************** */
/* void BVRShift( BitVect R, BitVect A, int n )                      */
/* Effectue : R = A >> n ;                                                */
/* ********************************************************************** */
void BVRShift (BitVect * R, BitVect * A, int n)
{
   int i;
   int WLmn;
   unsigned long temp;
   if ((R->n != A->n)) {
      printf ("Error in BVRShift(): Vectors of different sizes\n");
      exit (1);
   }

   for (i = 0; i < A->n; i++)
      R->vect[i] = A->vect[i];
   while (n >= 32) {
      for (i = A->n; i > 1; i--)
         R->vect[i - 1] = R->vect[i - 2];
      R->vect[0] = 0UL;
      n -= 32;
   }
   if (n > 0) {
      WLmn = WL - n;
      R->vect[A->n - 1] >>= n;
      for (i = A->n - 2; i >= 0; i--) {
         temp = R->vect[i] << WLmn;
         R->vect[i + 1] |= temp;
         R->vect[i] >>= n;
      }
   }
}


/* ********************************************************************** */
/* void BVLShiftSelf( BitVect R, int n )                             */
/* Effectue : R <<= n ;                                                   */
/* ********************************************************************** */
void BVLShiftSelf (BitVect * R, int n)
{
   int i;
   int WLmn;
   unsigned long temp;

   while (n >= 32) {
      for (i = 1; i < R->n; i++)
         R->vect[i - 1] = R->vect[i];
      R->vect[R->n - 1] = 0UL;
      n -= 32;
   }
   if (n > 0) {
      WLmn = WL - n;
      R->vect[0] <<= n;
      for (i = 1; i < R->n; i++) {
         temp = R->vect[i] >> WLmn;
         R->vect[i - 1] |= temp;
         R->vect[i] <<= n;
      }
   }
}


/* ********************************************************************** */
/* void BVLS1Self( BitVect R )                                            */
/* Effectue : R <<= 1 ;                                                   */
/* Version specialisee de BVLShiftSelf pour utilisation frequente.        */
/* ********************************************************************** */
void BVLS1Self (BitVect * R)
{
   int i;
   R->vect[0] <<= 1;
   for (i = 1; i < R->n; i++) {
      if (R->vect[i] & MC)
         R->vect[i - 1] |= 0x1UL;
      R->vect[i] <<= 1;
   }
}


/* ********************************************************************** */
/* void BVRShiftSelf( BitVect R, int n )                             */
/* Effectue : R >>= n ;                                                   */
/* ********************************************************************** */
void BVRShiftSelf (BitVect * R, int n)
{
   int i;
   int WLmn;
   unsigned long temp;

   while (n >= 32) {
      for (i = R->n - 1; i > 0; i--)
         R->vect[i] = R->vect[i - 1];
      R->vect[0] = 0UL;
      n -= 32;
   }
   if (n > 0) {
      WLmn = WL - n;
      R->vect[R->n - 1] >>= n;
      for (i = R->n - 2; i >= 0; i--) {
         temp = R->vect[i] << WLmn;
         R->vect[i + 1] |= temp;
         R->vect[i] >>= n;
      }
   }
}


/* ********************************************************************** */
/* void invertBV(Bitvect A)                                               */
/* fait A ~=A                                                             */
/* ********************************************************************** */
void InverseBV (BitVect * A)
{
   int i;
   for (i = 0; i < A->n; i++)
      A->vect[i] = ~A->vect[i];
}


/* ********************************************************************** */
/* lebool CheckCD( BitVect ds1, BitVect ds2 )                            */
/* Verifie si les ensembles ds1 et ds2 ont des bits communs.              */
/* ********************************************************************** */
lebool VerifBitsCommuns (BitVect * ds1, BitVect * ds2)
{
   int i;
   unsigned long temp = 0UL;
   if ((ds1->n != ds2->n)) {
      printf ("Error in VerifBitsCommuns(): Vectors of different sizes\n");
      exit (1);
   }
   for (i = 0; i < ds1->n; i++)
      temp |= (ds1->vect[i] & ds2->vect[i]);
   if (temp)
      return TRUE;
   else
      return FALSE;
}


/* -------------------------------------------- */
/* Fonctions pour la manipulation des matrices. */
/* -------------------------------------------- */
void BVCanonic (BitVect * A, int l)
{
   int n;
   PutBVToZero (A);
   n = l / WL;
   if (n > A->n) {
      printf
         ("Error in  BVCanonic(): vector A is not long enough to store  BVCanonic[%d].\n",
         l);
      exit (1);
   }
   A->vect[n] = 0x80000000UL >> (l - n * WL);
}

void Mask (BitVect * A, int l)
{

   InvMask (A, l);
   InverseBV (A);
}

void InvMask (BitVect * A, int l)
{
   AllOnes (A);
   BVRShiftSelf (A, l);
}

void AllOnes (BitVect * A)
{
   int i;
   for (i = 0; i < A->n; i++)
      A->vect[i] = 0xffffffffUL;
}

void AllocBV (BitVect * A, int l)
{
   int n;
   n = (l - 1) / WL + 1;
   A->vect = (unsigned long *) calloc ((size_t) n, sizeof (unsigned long));
   A->n = n;
}

void FreeBV (BitVect * A)
{
   if (A->vect != NULL) {
      free (A->vect);
   }
   A->n = 0;
}

void AllocMat (Matrix * m, int nblines, int l, int t)
{
   int i, j;
   m->lignes = (BitVect **) calloc ((size_t) nblines, sizeof (BitVect *));
   for (i = 0; i < nblines; i++) {
      if (!(m->lignes[i] = (BitVect *) calloc ((size_t) t, sizeof (BitVect)))) {
         printf ("\n*** Memoire insuffisante pour AllocMat() ! nl=%d***\n",
            nblines);
         exit (1);
      }
      for (j = 0; j < t; j++)
         AllocBV (&(m->lignes[i][j]), l);

   }
   m->nblignes = nblines;
   m->t = t;
   m->l = l;

}


/* ********************************************************************** */
/* void FreeSpace( Matrix m, int nl )                               */
/* Libere l'espace des nl vecteurs de la matrice m.                       */
/* ********************************************************************** */
void FreeMat (Matrix * m)
{
   int i, j;
   for (i = 0; i < m->nblignes; i++) {
      for (j = 0; j < m->t; j++)
         FreeBV (&(m->lignes[i][j]));
      free (m->lignes[i]);
   }
   free (m->lignes);
   m->nblignes = 0;
   m->l = 0;
   m->t = 0;

}


/* ********************************************************************** */
/* void CopyMat( Matrix m, Matrix ms, int nl, int t )                   */
/* Cette procedure sert a copier les t premiers BitVect des nl premieres  */
/* lignes de la matrice ms dans la matrice de travail m.                  */
/* ********************************************************************** */
void CopyMat (Matrix * m, Matrix * ms, int nl, int t)
{
   int i, j;
   if (m == NULL) {
      AllocMat (m, ms->nblignes, ms->l, ms->t);
   } else if ((ms->nblignes < nl) || (ms->t < t)) {
      printf ("Error in CopyMat(): source matrix too small %d\n",
         ms->nblignes / ms->t);
      exit (1);
   } else if ((m->nblignes < nl) || (m->t < t)) {
      printf ("Error in CopyMat(): destination matrix too small\n");
      exit (1);
   }
   for (i = 0; i < nl; i++)
      for (j = 0; j < t; j++) {
         CopyBV (&(m->lignes[i][j]), &(ms->lignes[i][j]));
      }
}


/* ********************************************************************** */
/* void CopyNTupleMat( Matrix m, Matrix ms, int nl,     */
/*                     int *colonnes, int t )                   */
/* Cette procedure sert a copier les t-1 BitVect indiqu'es par le vecteur  */
/* *colonnes plus la colonne 0 des nl premieres lignes de la matrice ms   */
/* dans la matrice de travail m.       */
/* ********************************************************************** */
void CopyNTupleMat (Matrix * m, Matrix * ms, int nl, int *colonnes, int t)
{

   int i, j, k, n;

   if (m == NULL)
      AllocMat (m, ms->nblignes, ms->l, t);
   else {
      if ((ms->nblignes != m->nblignes) || (ms->l != m->l))
         printf ("Error in CopieNTupleMat(): matrices of different sizes\n");
   }
   n = (ms->l - 1) / WL;
   for (i = 0; i < nl; i++) {
      for (k = 0; k <= n; k++)
         (m->lignes[i])[0].vect[k] = (ms->lignes[i])[0].vect[k];
      for (j = 1; j < t; j++)
         for (k = 0; k <= n; k++)
            (m->lignes[i])[j].vect[k] =
               (ms->lignes[i])[colonnes[j - 1]].vect[k];

   }
}


/* ********************************************************************** */
/* void SwapVect( Matrix m, int i, int j )                     */
/* Pour interchanger les lignes i et j de la matrice de travail m.        */
/* ********************************************************************** */
void ExchangeVect (Matrix * m, int i, int j)
{
   BitVect *p;
   if (i != j) {
      p = m->lignes[i];
      m->lignes[i] = m->lignes[j];
      m->lignes[j] = p;
   }
}


void TransposeMatrices (Matrix * T, Matrix * M, int mmax, int smax, int L)
{

   int s, l, m;

   for (s = 0; s < smax; s++)
      for (l = 0; l < L; l++) {
         PutBVToZero (&T->lignes[l][s]);
         for (m = 0; m < mmax; m++) {
            /* printf("m=%d l=%d s=%d\n",m,l,s);fflush(stdout); */
            if (M->lignes[m][s].vect[0] & (0x80000000UL >> l)) {
               T->lignes[l][s].vect[0] |= (0x80000000UL >> m);
            }
         }
      }
}


/* ********************************************************************** */
/* void XorVect( Matrix m,                                               */
/*               int r, int s, int min, int max )     */
/* Effectue un Xor entre la s-eme et la r-eme ligne de la matrice de      */
/* travail m pour les colonnes (BitVect) min a max-1 seulement.           */
/* Le resultat est mis dans la r-eme ligne. ( m[r] ^= m[s] )              */
/* ********************************************************************** */
void XorVect (Matrix * m, int r, int s, int min, int max)
{
   int j;
   for (j = min; j < max; j++)
      XORBVSelf (&(m->lignes[r][j]), &(m->lignes[s][j]));
}


/* ********************************************************************** */
/* void displaymat(Matrix m, int t, int l, int kg)        */
/* Affiche la matrice m sur kg lignes par t x l colonnes  		  */
/* ********************************************************************** */
void DispMat (Matrix * m, int t, int l, int kg, lebool mathematica)
{
   int i, j;

   i = kg;

   printf ("\n");
   if (mathematica)
      printf ("{");
   for (i = 0; i < kg; i++) {
      if (!mathematica)
         printf ("[");
      for (j = 0; j < t; j++) {
         DispBitVect (&(m->lignes[i][j]), l, mathematica);
      }
      if (mathematica) {
         if (i != kg - 1)
            printf (",\n");
         else
            printf ("}\n");
      } else
         printf ("]\n");
   }
   printf ("\n\n");

}


/* ********************************************************************** */
/* void displaybitvect(BitVect A, int l)    			  */
/* Affiche le BitVect A sur l bits seulement    			  */
/* ********************************************************************** */
void DispBitVect (BitVect * A, int l, int mathematica)
{
   int j;
   unsigned Un;
   Un = 1UL;
   j = 0;
   if (mathematica) {
      printf ("{");
      while (j < l - 1) {
         printf ("%ld,",
            (A->vect[j / WL] >> (((WL * A->n) - j - 1) % WL)) & Un);
         j++;
      }
      printf ("%ld}", (A->vect[j / WL] >> (((WL * A->n) - j - 1) % WL)) & Un);
   } else
      while (j < l) {
         printf ("%ld",
            (A->vect[j / WL] >> (((WL * A->n) - j - 1) % WL)) & Un);
         j++;
      }
}


/* ********************************************************************** */
/* MultMatrixByMatrix ( Matrix *A, Matrix *B, Matrix *C)        	  */
/* Fait la multiplication matricielle A = B x C    			  */
/* ********************************************************************** */

void MultMatrixByMatrix (Matrix * A, Matrix * B, Matrix * C)
{
   int i, j;
   if (B->l != C->nblignes) {
      printf ("Tailles de matrices non-compatibles, kaput.\n");
      exit (1);
   }

   if (A->nblignes != B->nblignes || A->l != C->l) {
      printf ("Matrice preallouee de mauvaise taille.\n");
      exit (1);
   }

   for (i = 0; i < A->nblignes; i++)
      PutBVToZero (A->lignes[i]);

   for (i = 0; i < B->nblignes; i++)
      for (j = 0; j < B->l; j++) {
         if (ValBitBV (B->lignes[i], j))
            XORBVSelf (A->lignes[i], C->lignes[j]);
      }
}


/* ********************************************************************** */
/* MatrixTwoPow ( Matrix *A, Matrix *B, unsigned int e)          	  */
/* Fait l'exponentiation : A = B^(2^e)          			  */
/* ********************************************************************** */

void MatrixTwoPow (Matrix * A, Matrix * B, unsigned int e)
{
   unsigned int i;
   Matrix tempMatrix;
   Matrix *AA = &tempMatrix;

   if (B->nblignes != B->l) {
      printf ("Matrice non carree.\n");
      exit (1);
   }
   if (A->nblignes != B->nblignes || A->l != B->l) {
      printf ("Matrice preallouee de mauvaise taille.\n");
      exit (1);
   }

   AllocMat (AA, B->nblignes, B->l, 1);

   if (e == 0) {
      CopyMat (A, B, B->nblignes, 1);
      return;
   }
   /* A = B^(2^1) */
   MultMatrixByMatrix (A, B, B);

   for (i = 1; i < e - 1; i += 2) {
     /* AA = A * A */
      MultMatrixByMatrix (AA, A, A);

      /* A = AA * AA */
      MultMatrixByMatrix (A, AA, AA);
   }

   if (i == e - 1) {
     /* AA = A * A */
      MultMatrixByMatrix (AA, A, A);
      CopyMat (A, AA, AA->nblignes, 1);
   }

   FreeMat (AA);
}

/* ********************************************************************** */
/* MatrixPow ( Matrix *A, Matrix *B, unsigned int e)            	  */
/* Fait l'exponentiation : A = B^e               			  */
/* ********************************************************************** */

#ifdef USE_LONGLONG
void MatrixPow (Matrix * A, Matrix * B, longlong e)
#else
void MatrixPow (Matrix * A, Matrix * B, long e)
#endif
{
   int i;
   Matrix C;
   Matrix D;

   if (B->nblignes != B->l) {
      printf ("Matrice non carree.\n");
      exit (1);
   }

   if (A->nblignes != B->nblignes || A->l != B->l) {
      printf ("Matrice preallouee de mauvaise taille.\n");
      exit (1);
   }

   AllocMat (&C, B->nblignes, B->l, 1);

   if (e < 0) {
      InverseMatrix (&C, B);
      MatrixPow (A, &C, -e);
      FreeMat (&C);
      return;
   }
   AllocMat (&D, B->nblignes, B->l, 1);

   /* A = I */
   for (i = 0; i < A->nblignes; i++)
      BVCanonic (A->lignes[i], i);

   /* C = B^1 */
   CopyMat (&C, B, B->nblignes, 1);

   while (e) {
      if (e & 1) {
         CopyMat (&D, A, B->nblignes, 1);
         MultMatrixByMatrix (A, &D, &C);
      }

      e >>= 1;

      if (e) {
         CopyMat (&D, &C, B->nblignes, 1);
         MultMatrixByMatrix (&C, &D, &D);
      }
   }

   FreeMat (&C);
   FreeMat (&D);
}


/* FIN vectorsF2.c */
