/*****************************************************************************/
/*             A Class to Manipulate Hermite Polynomials                      */
/*****************************************************************************/

#include <cmath>
#include <cstdlib>
#include "VECMAT.h"
#include "Teuchos_Array.hpp"

#ifndef CHAOS_H
#define CHAOS_H

int nterm(int x, int y);        /* function to determine the number of terms corresponding to each specific order */

class Chaos {

 private:
  int ndim;                    /* number of dimensions in the Chaos Expansionn */
  int order;                   /* order of the expansion */

 public:
  Chaos(int dim, int order);
  int tnterms();   /* total number of terms in a chaos expansion */
  void HermiteToPowers(int &maxorder, George::matrix &HP); /* generates the different orders of a 1d hermite polynomial */
  void ScHermiteToPowers(int &maxorder, double &scale, George::matrix &HP);
  void mult1dHermite(George::vector &A, George::vector &B, George::vector &C);
  void OrthoOrder(George::matrix &D);
};

inline Chaos::Chaos(int nd, int od) {
  ndim = nd;
  order = od;
}
inline int Chaos::tnterms(){
  int x = 1;
  for (int i=1; i<=order; i++)
    x += nterm(i, ndim);
  return x;
}

inline void Chaos::HermiteToPowers(int &mr, George::matrix &icoef){

  const int maxorder = mr;

  for (int i=0; i<=maxorder; i++)
    {
      for (int j=0; j<=maxorder; j++)
        {
          icoef[i][j] = 0;
        }
    }
  icoef [1][1] = 1;
  icoef [2][2] = 1;

  for (int k=2; k<=maxorder; k++)
    {
      for (int i=1; i<=maxorder; i++)
        {icoef[k][i] = icoef[k-1][i-1];
        }
      for (int i=0; i<=maxorder; i++)
        {icoef[k][i] = icoef[k][i]-icoef[k-2][i]*(k-2);
        }
    }
}

inline void Chaos::ScHermiteToPowers(int &mr,double & sc, George::matrix &icoef){

  const int maxorder = mr;
  //double scale = sc;

  for (int i=0; i<=maxorder; i++)
    {
      for (int j=0; j<=maxorder; j++)
        {
          icoef[i][j] = 0;
        }
    }
  icoef [1][1] = 1;
  icoef [2][2] = 1;

  for (int k=2; k<=maxorder; k++)
    {
      for (int i=1; i<=maxorder; i++)
        {icoef[k][i] = icoef[k-1][i-1];
        }
      for (int i=0; i<=maxorder; i++)
        {icoef[k][i] = icoef[k][i]-icoef[k-2][i]*(k-2);
        }
    }
  for(int i=1; i<=maxorder; i++)
    {
      for(int j=1; j<=maxorder; j++)
        {
          if(j<i)
            icoef[i][j] *= pow(sc,(i-j));
        }
    }
}

inline void Chaos::mult1dHermite(George::vector &A, George::vector &B, George::vector &C){
  int noi = A.getsize();
  int noj = B.getsize();
  int nok = C.getsize();

  // arrays containing the coef. of each of the polynomials

  for(int i=0; i<nok; i++)
    C[i] = 0;

  for (int i=0; i<noi; i++)
    {
      for(int j=0; j<noj; j++)
        {
          C[i+j] += A[i]*B[j];
        }
    }
}
inline void Chaos::OrthoOrder(George::matrix &D){

  int nterms;
  int  maxterms = nterm(order, ndim);

  /* KL - changed from C to Teuchos arrays to be Ansi C++ compliant */
  /* KL 13 June 2008 -- removed two-argument ctors because of a compilation
   * problem under gcc 4.3.0. */
  //  Teuchos::Array<Teuchos::Array<int> > c(maxterms+1, order+1);
  // Teuchos::Array<Teuchos::Array<int> > c1(maxterms+1, ndim+1);
  Teuchos::Array<Teuchos::Array<int> > c(maxterms+1);
  for (int i=0; i<=maxterms; i++) c[i].resize(order+1);

  Teuchos::Array<Teuchos::Array<int> > c1(maxterms+1);
  for (int i=0; i<=maxterms; i++) c1[i].resize(ndim+1);



  int count = 0;
  int count1 = 0;

  for (int iorder=1; iorder<=order; iorder++)
    {
      nterms = nterm(iorder, ndim);

      for (int i=1; i<=iorder; i++)
        {
          for(int j=1; j<=nterms; j++)
            {
              c[j][i] = 0;
            }
        }

      int iterm = 1;
      count++;

      for(int i=1; i<=iorder; i++)
        c[iterm][i] = 1;


      for(int i=1; i<=ndim; i++)
        {
          c1[iterm][i] = 0;
          for(int j=1; j<=iorder; j++)
            {
              if (c[iterm][j]== i)
                c1[iterm][i]+=1;
            }
        }


      while (iterm < nterms)
        {
          iterm++;
          int iswitch = 0;
          count++;

          for(int i=iorder; i>=1; i--)
            {
              if(((c[iterm-1][i]+1)<=ndim)&&(iswitch==0))
                {
                  iswitch = 1;
                  c[iterm][i] = c[iterm-1][i] + 1;

                  if(i < iorder)
                    {
                      for(int j=i+1; j<=iorder; j++)
                        c[iterm][j] = c[iterm][i];
                    }

                  if (i>1)
                    {
                      for(int j=1; j<=i-1; j++)
                        c[iterm][j] = c[iterm-1][j];
                    }
                }
            }



          for (int i=1;  i<=ndim; i++)
            {
              c1[iterm][i] = 0;
              for(int j=1; j<=iorder; j++)
                {
                  if( c[iterm][j] == i)
                    c1[iterm][i]+=1;
                }
            }
        }


      for (int it = 1; it<=nterms; it++)
        {
          count1++;
          for(int j=1; j<=ndim; j++)
            {
              D[count1][j] = c1[it][j];
            }
        }
    }
}

/* defining the function nterm */

inline int nterm(int x,int y)
{

  int fac = 1;
  int nterm = 1;

  for (int i=0; i<x; i++)
    {
      if (i>0)
        fac = fac * i;
      nterm = nterm*(y+i);
    }
  fac = fac*x;
  if (fac>0)
    nterm = nterm/fac;
  else
    nterm = 1;

  return nterm;
}

#endif
