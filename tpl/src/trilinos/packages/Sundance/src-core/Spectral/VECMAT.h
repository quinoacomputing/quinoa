// Classes Vector and Matrix

#ifndef VECMAT
#define VECMAT

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>


namespace George
{

  class matrix;                                       // Incomplete Declaration

  class vector{
    friend class matrix;
  private:
    int size;
    double *vec;
    int range(int);                                        // index check
  public:
    vector(int);                                           // init to zero elems
    vector(const vector&);                                 // init to given vec
    vector(const double *, int);                      // init to given array
    ~vector() {delete [] vec;}
    double& operator[](int i) {return vec[range(i)];}
    vector& operator+()        {return *this;}          // unary plus
    vector& operator=(const vector&);
    vector& operator+=(const vector&);
    vector& operator-=(const vector&);
    vector& operator*=(double);                         // multiply by double
    vector& operator/=(double);                         // divide by double
    int getsize() {return size;}
    void swap(int, int);
    friend vector operator-(const vector&);
    friend vector operator+(const vector&, const vector&);
    friend vector operator-(const vector&, const vector&);
    friend vector operator*(const vector&, double);
    friend vector operator*(double, const vector&);
    friend vector operator/(const vector&, double);
    friend double scalar(const vector&, const vector&);    // scalar product
    friend double norm(const vector&);                    // Euclidean norm
    friend double norminf(const vector&);                 // infinite norm
    friend vector operator*(const matrix&, const vector&); // v = m.u
    friend vector operator*(const vector&, const matrix&); // v = u.m
    friend matrix operator*(const vector&, const vector&); // m = u.v
    friend matrix operator*(const matrix&, const matrix&); // m = m1.m2
    friend double norm(const matrix&);                    // Euclidean norm
    friend double norminf(const matrix);                 // infinite norm

  };
  class matrix{
  private:
    int numrows;
    int numcols;
    vector **mat;
    int range(int);                                     //row index check

  public:
    matrix(int, int);                                  // rectangular matrix
    matrix(int);                                       //square matrix
    matrix(const matrix&);                            // init to given matrix
    ~matrix();
    vector& operator[](int i) {return *mat[range(i)];}
    matrix& operator+() {return *this;}
    matrix& operator=(const matrix&);
    matrix& operator+=(const matrix&);
    matrix& operator-=(const matrix&);
    matrix& operator*=(double);
    matrix& operator/=(double);
    int getnumrows() {return numrows;}
    int getnumcols() {return numcols;}
    int getsize();                               // square matrix only
    void swap(int, int);
    matrix transpose();
    friend matrix operator-(const matrix&);
    friend matrix operator+(const matrix&, const matrix&);
    friend matrix operator-(const matrix&, const matrix&);
    friend matrix operator*(const matrix&, const matrix&);  // m = m1.m2
    friend matrix operator*(const matrix&, double);        // m = m1.d
    friend matrix operator*(double,const matrix&);        // m = d.m1
    friend matrix operator/(const matrix&, double);       // m = m1/d
    friend vector operator*(const vector&, const matrix&); // v= u.m
    friend vector operator*(const matrix&, const vector&); //v = m.u
    friend matrix operator*(const vector&, const vector&);  // m = u.v
    friend double norm(const matrix&);                  // Euclidean norm
  };

  // inline functions

  inline int vector::range(int i)
    {if (i<0 || i>=size)
        printf("ERROR!! vector index out of range!\n");
      return i;
    }

  inline int matrix::range(int i)
    {if (i<0 || i>=numrows)
        printf("ERROR matrix row index out of range!\n");
      return i;
    }
  inline vector operator*(double d, const vector &v)
    { return v*d;}
  inline matrix operator*(double d, const matrix &m)
    { return m*d; }


  inline vector::vector(int n)
    {size = n;
      vec = new double [size];
      if(!vec)
        printf("ERROR!! allocation failure in vector::vector(int)!\n");
      for (int i=0; i<size; i++)
        vec[i] = 0;
    }
  inline vector::vector(const double *a, int n)
    {size = n;
      vec = new double [size];
      if(!vec)
        printf("ERROR!! allocation failure in vector::vector(double*, int)!\n");
      for (int i=0; i<size; i++)
        vec[i] = a[i];
    }
  inline vector::vector(const vector &v)
    {size = v.size;
      vec = new double[size];
      if(!vec)
        printf("ERROR allocation failure in vector::vector(vector&)!\n");
      for (int i=0; i<size; i++)
        vec[i] = v.vec[i];
    }
  inline vector& vector::operator=(const vector &v)
    {if (size != v.size)
        printf("ERROR !! diff size in vector& vector::op=(const vector&)!\n");
      for (int i=0; i<size; i++)
        vec[i] = v.vec[i];
      return *this;
    }
  inline vector& vector::operator+=(const vector &v)
    {if (size != v.size)
        printf("ERROR diff size in vector& vector::op=(const vector&)!\n");
      for (int i=0; i<size; i++)
        vec[i] += v.vec[i];
      return *this;
    }
  inline vector& vector::operator-=(const vector &v)
    {if (size != v.size)
        printf("ERROR diff size in vector& vector::op=(const vector&)!\n");
      for (int i=0; i<size; i++)
        vec[i] -= v.vec[i];
      return *this;
    }
  inline vector& vector::operator*=(double x)
    {for (int i=0; i<size; i++)
        vec[i] *= x;
      return *this;
    }
  inline vector& vector::operator/=(double x)
    {for (int i=0; i<size; i++)
        vec[i] /= x;
      return *this;
    }
  inline void vector::swap(int i, int j)
    {double tmp = vec[range(i)];
      vec[i] = vec[range(j)];
      vec[j] = tmp;
    }
  inline vector operator-(const vector &v)
    {int n = v.size;
      vector u(n);
      for (int i=0; i<n; i++)
        u.vec[i]=-v.vec[i];
      return u;
    }
  inline vector operator+(const vector &v1, const vector &v2)
    {int n = v1.size;
      if(v1.size != v2.size)
        printf("error different vector sizes\n");
      vector v(n);
      for (int i=0; i<n; i++)
        v.vec[i] = v1.vec[i] + v2.vec[i];
      return v;
    }
  inline vector operator-(const vector &v1, const vector &v2)
    {int n = v1.size;
      if(v1.size != v2.size)
        printf("error different vector sizes\n");
      vector v(n);
      for (int i=0; i<n; i++)
        v.vec[i] = v1.vec[i] - v2.vec[i];
      return v;
    }
  inline vector operator*(const vector &v, double d)
    {vector vd = v;
      vd *= d;
      return vd;
    }
  inline vector operator/(const vector &v, double d)
    {vector vd = v;
      vd /= d;
      return vd;
    }
  inline double scalar(const vector &u, const vector &v)
    {double t=0;
      int n= u.size;
      if(u.size!= v.size)
        printf("error different vector sizes\n");
      for (int i=0; i<n; i++)
        t += u.vec[i]*v.vec[i];
      return t;
    }
  inline double norm(const vector &v)
    {int n = v.size;
      const double eps =  1.0E-15;
      double t = fabs(v.vec[0]);
      for (int i=0; i<n; i++)
        { double vi = fabs(v.vec[i]);
          if (t<eps) t=vi;
          else
            {double x = vi/t;
              t *= sqrt(1+x*x);
            }
        }
      return t;
    }
  inline double norminf(const vector &v)
    {int n = v.size;
      double t = 0.0;
      for (int i=0; i<n; i++)
        {double vi = fabs(v.vec[i]);
          if (vi > t ) t= vi;
        }
      return t;
    }

  inline matrix::matrix(int nrows, int ncols)
    { numrows = nrows;
      numcols = ncols;
      mat = new vector* [numrows];
      if (!mat)
        printf("error row allocation failure\n");
      for (int i = 0; i<numrows; i++)
        {mat[i] = new vector(numcols);
          if(!mat[i])
            printf("error column allocation failure\n");
        }
    }
  inline matrix::matrix(int n)
    { numrows = n;
      numcols = n;
      //std::cerr << "n=" << n << std::endl;
      mat = new vector* [numrows];
      if (!mat)
        printf("error row allocation failure\n");
      for (int i = 0; i<numrows; i++)
        {mat[i] = new vector(numcols);
          if(!mat[i])
            printf("error column allocation failure\n");
        }
      //std::cerr << "done matrix ctor" << std::endl;
    }
  inline matrix::matrix(const matrix &m)
    {numrows = m.numrows;
      numcols = m.numcols;
      mat = new vector* [numrows];
      if (!mat)
        printf("error row allocation failure\n");
      for (int i = 0; i<numrows; i++)
        {mat[i] = new vector(numcols);
          if(!mat[i])
            printf("error column allocation failure\n");
        }
      for (int i=0; i<numrows; i++)
        *mat[i] = *m.mat[i];
    }
  inline matrix::~matrix()
    {for (int i=numrows; i>0; --i)
        delete mat[i-1];
      delete [] mat;
    }
  inline matrix& matrix::operator=(const matrix &m)
    {if (m.numrows != numrows || numcols != numcols)
        printf("error diff sizes\n");
      for(int i= 0; i< numrows; i++)
        *mat[i] = *m.mat[i];
      return *this;
    }
  inline matrix& matrix::operator+=(const matrix &m)
    {if (m.numrows != numrows || numcols != numcols)
        printf("error diff sizes\n");
      for(int i= 0; i< numrows; i++)
        *mat[i] += *m.mat[i];
      return *this;
    }
  inline matrix& matrix::operator-=(const matrix &m)
    {if (m.numrows != numrows || numcols != numcols)
        printf("error diff sizes\n");
      for(int i= 0; i< numrows; i++)
        *mat[i] -= *m.mat[i];
      return *this;
    }
  inline matrix& matrix::operator*=(double x)
    {for (int i=0; i<numrows; i++)
        {for (int j=0; j<numcols; j++)
            mat[i]->vec[j] *=x;
        }
      return *this;
    }


  inline matrix& matrix::operator/=(double x)
    {for (int i=0; i<numrows; i++)
        {for (int j=0; j<numcols; j++)
            mat[i]->vec[j] /=x;
        }
      return *this;
    }

  inline int matrix::getsize()
    {if (numrows != numcols)
        printf("error getsize() requires a square matrix\n");
      return numrows;
    }
  inline void matrix::swap(int i, int j)
    {vector *tmp = mat[range(i)];
      mat[i] = mat [ range(j)];
      mat[j] = tmp;
    }
  inline matrix matrix::transpose()
    { int p = numrows;
      int q = numcols;
      matrix mt(q,p);
      for (int i = 0; i<q; i++)
        {for (int j=0; j<p; j++)
            mt.mat[i]->vec[j] = mat[j] ->vec[i];
        }
      return mt;
    }
  inline matrix operator-(const matrix &m)
    {int p = m.numrows;
      int q = m.numcols;
      matrix mt(q,p);
      for (int i=0; i<p; i++)
        {*mt.mat[i] = - *m.mat[i];
        }
      return mt;
    }
  inline matrix operator+(const matrix &m1, const matrix &m2)
    {if (m1.numrows != m2.numrows || m1.numcols != m2.numcols)
        printf("error diff sizes\n");
      matrix mt = m1;
      mt += m2;
      return mt ;
    }
  inline matrix operator-(const matrix &m1, const matrix &m2)
    {if (m1.numrows != m2.numrows || m1.numcols != m2.numcols)
        printf("error diff sizes\n");
      matrix mt = m1;
      mt -= m2;
      return mt ;
    }
  inline matrix operator*(const matrix &m, double d)
    {matrix mt = m;
      mt *=d;
      return mt;
    }
  inline matrix operator/(const matrix &m, double d)
    {matrix mt = m;
      mt /=d;
      return mt;
    }
  inline double norm(const matrix &m)
    {double t = 0;
      int nr = m.numrows;
      int nc = m.numcols;
      for (int i=0; i<nr ; i++)
        {for (int j=0; j<nc; j++)
            {double mij = m.mat[i]->vec[j];
              t+= mij*mij;
            }
        }
      return sqrt(t);
    }

  inline vector operator*(const matrix &m, const vector &v)
    {int nr = m.numrows;
      if (m.numcols!= v.size)
        printf("error diff sizes\n");
      vector u(nr);
      for (int i=0; i<nr; i++)
        u[i] = scalar(*m.mat[i],v);
      return u;
    }
  inline vector operator*(const vector &v, const matrix &m)
    {int nr = m.numrows;
      int nc = m.numcols;
      if (v.size != nr)
        printf("error diff sizes\n");
      vector u(nc);
      for (int i =0; i<nc; i++)
        {double t = 0;
          for (int j=0; j<nr; j++)
            t+= v.vec[j]*m.mat[j]->vec[i];
          u.vec[i] = t;
        }
      return u;
    }
  inline matrix operator*(const vector &u, const vector &v)
    {int nr = u.size;
      int nc = v.size;
      matrix m(nr, nc);
      for (int i=0; i<nr; i++)
        {for (int j=0; j<nc; j++)
            m.mat[i]->vec[j] = u.vec[i] * v.vec[j];
        }
      return m;
    }
  inline matrix operator*(const matrix &m1, const matrix &m2)
    {int nr1 = m1.numrows;
      int nc1 = m1.numcols;
      int nr2 = m2.numrows;
      int nc2 = m2.numcols;
      if (nc1 != nr2 )
        printf("error matrices are of different sizes\n");
      matrix m(nr1,nc2);
      for (int i=0; i<nr1; i++)
        {for(int k=0; k<nc2; k++)
            {double t=0;
              for(int j=0; j<nc1; j++)
                t+= m1.mat[i]->vec[j]*m2.mat[j]->vec[k];
              m.mat[i]->vec[k] = t;
            }
        }
      return m;
    }
}

#endif
