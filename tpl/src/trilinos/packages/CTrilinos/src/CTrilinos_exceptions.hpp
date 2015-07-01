
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


/*! @file CTrilinos_exceptions.hpp
 * @brief Defines exceptions thrown by CTrilinos. */


#ifndef CTRILINOS_EXCEPTIONS_HPP
#define CTRILINOS_EXCEPTIONS_HPP


#include "CTrilinos_config.h"


#include "Teuchos_Exceptions.hpp"


namespace CTrilinos {


/*! exception indicating wrong object type encountered */
class CTrilinosTypeMismatchError : public Teuchos::ExceptionBase
{
  public:
    CTrilinosTypeMismatchError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};

/*! exception indicating wrong object type encountered */
class CTrilinosInvalidTypeError : public Teuchos::ExceptionBase
{
  public:
    CTrilinosInvalidTypeError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};

/*! exception indicating invalid cast attempted */
class CTrilinosConstCastError : public Teuchos::ExceptionBase
{
  public:
    CTrilinosConstCastError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};

/*! exception indicating wrong object table accessed */
class CTrilinosWrongTableError : public Teuchos::ExceptionBase
{
  public:
    CTrilinosWrongTableError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};

/*! exception indicating some other error */
class CTrilinosMiscException : public Teuchos::ExceptionBase
{
  public:
    CTrilinosMiscException(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};


} // namespace CTrilinos


#endif

