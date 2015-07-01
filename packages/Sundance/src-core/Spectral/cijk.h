/*****************************************************************************/
/*          A Class to to Calculate the Expectated Value                     */
/*              of a product of Hermite polynomials                          */
/*****************************************************************************/

#include <cmath>
#include <cstdlib>
#include "VECMAT.h"
#include "chaos.h"
#include "Teuchos_Array.hpp"

#ifndef CIJK_H
#define CIJK_H


class cijk {

private:
  int ndim;         /* number of dimensions in expansion */
  int order;        /* order of expansion */
  int maxorder, maxnterms;
  int i, j, k,l,m;
  void PC1(int mn, int n);
  double Orthorder[2500][100];
  double icoef[100][100];

public:
  cijk(int ndim, int order);
  double expectation(int i, int j, int k);
  double expectation4(int i, int j, int k, int l);
  double expectation5(int i, int j, int k, int l, int m);
};
inline void cijk::PC1(int mn, int n)
{
  Chaos PC1(mn,n);
  maxnterms = PC1.tnterms();
  maxorder = 11;
  George::matrix  Orth(maxnterms, ndim+1);
  George::matrix coef(maxorder+1);

  PC1.HermiteToPowers(maxorder, coef);
  PC1.OrthoOrder(Orth);

  for(int a=0; a<maxnterms; a++)
    for(int b=0; b<ndim+1; b++)
      Orthorder[a][b] = Orth[a][b];

  for(int a=0; a<maxorder+1; a++)
    for(int b=0; b<maxorder+1; b++)
      icoef[a][b] = coef[a][b];

}

inline cijk::cijk(int nd, int od){
  ndim = nd;
  order = od;
  PC1(ndim, order);
}
inline double cijk::expectation(int i, int j, int k){

  Chaos PC(ndim, order);
  maxnterms = PC.tnterms();

  if (i>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (j>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (k>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");

  George::vector A(ndim+1);
  George::vector B(ndim+1);
  George::vector C(ndim+1);

  for(int id=1; id<=ndim; id++)
  {
    A[id] = int(Orthorder[i][id]);
    B[id] = int(Orthorder[j][id]);
    C[id] = int(Orthorder[k][id]);
    if (i==0)
      A[id] = 0;
    if (j==0)
      B[id] =0;
    if (k==0)
      C[id] = 0;
  }

  double prod = 1;

  Teuchos::Array<double> Mu(ndim+1);
  Teuchos::Array<double> Sigma(ndim+1);

  for(int in=1; in<=ndim; in++)
  {
    Mu[in] = 0;
    Sigma[in] = 1;
  }

  for(int dm=1; dm<=ndim; dm++)
  {
    double m = Mu[dm];
    double s = Sigma[dm];

    double Momt[40];
    Momt[0] = 1;
    Momt[1] = m;
    Momt[2] = pow(m,2) + pow(s,2);
    for(int mg = 3; mg<40; mg++)
      Momt[mg] = m*Momt[mg-1]+((mg-1)*pow(s,2)*Momt[mg-2]);

    George::vector A1(maxorder);
    George::vector B1(maxorder);
    George::vector C1(maxorder);

    George::vector P1(2*maxorder);
    George::vector P2(3*maxorder);

    for(int in=1; in<=maxorder; in++)
    {
      A1[in-1] = (icoef[int(A[dm]+1)][in]);
      B1[in-1] = (icoef[int(B[dm]+1)][in]);
      C1[in-1] = (icoef[int(C[dm]+1)][in]);
    }

    PC.mult1dHermite(A1,B1,P1);
    PC.mult1dHermite(P1,C1,P2);

    double p1=0;
    for(int ig=0; ig<3*maxorder; ig++)
      p1 += P2[ig]*Momt[ig];

    prod *= p1;
  }

  return prod;

}

inline double cijk::expectation4(int i, int j, int k, int l){

  Chaos PC(ndim, order);
  maxnterms = PC.tnterms();

  if (i>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (j>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (k>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (l>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");

  George::vector A(ndim+1);
  George::vector B(ndim+1);
  George::vector C(ndim+1);
  George::vector D(ndim+1);

  for(int id=1; id<=ndim; id++)
  {
    A[id] = int(Orthorder[i][id]);
    B[id] = int(Orthorder[j][id]);
    C[id] = int(Orthorder[k][id]);
    D[id] = int(Orthorder[l][id]);

    if (i==0)
      A[id] = 0;
    if (j==0)
      B[id] =0;
    if (k==0)
      C[id] = 0;
    if(l==0)
      D[id] = 0;
  }

  double prod = 1;

  Teuchos::Array<double> Mu(ndim+1);
  Teuchos::Array<double> Sigma(ndim+1);

  for(int in=1; in<=ndim; in++)
  {
    Mu[in] = 0;
    Sigma[in] = 1;
  }

  for(int dm=1; dm<=ndim; dm++)
  {
    double m = Mu[dm];
    double s = Sigma[dm];

    double Momt[60];
    Momt[0] = 1;
    Momt[1] = m;
    Momt[2] = pow(m,2) + pow(s,2);
    for(int mg = 3; mg<60; mg++)
      Momt[mg] = m*Momt[mg-1]+((mg-1)*pow(s,2)*Momt[mg-2]);

    George::vector A1(maxorder);
    George::vector B1(maxorder);
    George::vector C1(maxorder);
    George::vector D1(maxorder);

    George::vector P1(2*maxorder);
    George::vector P2(3*maxorder);
    George::vector P3(4*maxorder);

    for(int in=1; in<=maxorder; in++)
    {
      A1[in-1] = (icoef[int(A[dm]+1)][in]);
      B1[in-1] = (icoef[int(B[dm]+1)][in]);
      C1[in-1] = (icoef[int(C[dm]+1)][in]);
      D1[in-1] = (icoef[int(D[dm]+1)][in]);
    }

    PC.mult1dHermite(A1,B1,P1);
    PC.mult1dHermite(P1,C1,P2);
    PC.mult1dHermite(P2,D1,P3);

    double p1=0;
    for(int ig=0; ig<4*maxorder; ig++)
      p1 += P3[ig]*Momt[ig];

    prod *= p1;
  }

  return prod;

}
inline double cijk::expectation5(int i, int j, int k, int l,int m){

  Chaos PC(ndim, order);
  maxnterms = PC.tnterms();

  if (i>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (j>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (k>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (l>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (m>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");

  George::vector A(ndim+1);
  George::vector B(ndim+1);
  George::vector C(ndim+1);
  George::vector D(ndim+1);
  George::vector E(ndim+1);

  for(int id=1; id<=ndim; id++)
  {
    A[id] = int(Orthorder[i][id]);
    B[id] = int(Orthorder[j][id]);
    C[id] = int(Orthorder[k][id]);
    D[id] = int(Orthorder[l][id]);
    E[id] = int(Orthorder[m][id]);

    if (i==0)
      A[id] = 0;
    if (j==0)
      B[id] =0;
    if (k==0)
      C[id] = 0;
    if(l==0)
      D[id] = 0;
    if(m==0)
      E[id] = 0;
  }

  double prod = 1;

  Teuchos::Array<double> Mu(ndim+1);
  Teuchos::Array<double> Sigma(ndim+1);

  for(int in=1; in<=ndim; in++)
  {
    Mu[in] = 0;
    Sigma[in] = 1;
  }

  for(int dm=1; dm<=ndim; dm++)
  {
    double m = Mu[dm];
    double s = Sigma[dm];

    double Momt[70];
    Momt[0] = 1;
    Momt[1] = m;
    Momt[2] = pow(m,2) + pow(s,2);
    for(int mg = 3; mg<70; mg++)
      Momt[mg] = m*Momt[mg-1]+((mg-1)*pow(s,2)*Momt[mg-2]);

    George::vector A1(maxorder);
    George::vector B1(maxorder);
    George::vector C1(maxorder);
    George::vector D1(maxorder);
    George::vector E1(maxorder);

    George::vector P1(2*maxorder);
    George::vector P2(3*maxorder);
    George::vector P3(4*maxorder);
    George::vector P4(5*maxorder);

    for(int in=1; in<=maxorder; in++)
    {
      A1[in-1] = (icoef[int(A[dm]+1)][in]);
      B1[in-1] = (icoef[int(B[dm]+1)][in]);
      C1[in-1] = (icoef[int(C[dm]+1)][in]);
      D1[in-1] = (icoef[int(D[dm]+1)][in]);
      E1[in-1] = (icoef[int(E[dm]+1)][in]);
    }

    PC.mult1dHermite(A1,B1,P1);
    PC.mult1dHermite(P1,C1,P2);
    PC.mult1dHermite(P2,D1,P3);
    PC.mult1dHermite(P3,E1,P4);


    double p1=0;
    for(int ig=0; ig<5*maxorder; ig++)
      p1 += P4[ig]*Momt[ig];

    prod *= p1;
  }

  return prod;

}

#endif
