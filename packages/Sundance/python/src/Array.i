// -*- c++ -*-

%{
  // System includes
#include <Python.h>
#include "Teuchos_Array.hpp"
  %}


namespace Teuchos
{

template <typename T> class Array
{
public:
  Array(int n);

  int size() const ;

  void resize(int newsize);

  %extend
   {
     T __getitem__(int i) 
     {
       return self->operator[](i);
     }

     void __setitem__(int i, const T &x) 
     {
       self->operator[](i) = x;
     }

     std::string __str__()
     {
       return self->toString();
     }
   }
};




}

