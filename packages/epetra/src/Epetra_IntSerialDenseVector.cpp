
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#include "Epetra_IntSerialDenseVector.h"

//=============================================================================
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector()
  : Epetra_IntSerialDenseMatrix()
{
	SetLabel("Epetra::IntSerialDenseVector");
}

//=============================================================================
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(int length)
  : Epetra_IntSerialDenseMatrix(length, 1)
{
	SetLabel("Epetra::IntSerialDenseVector");
}

//=============================================================================
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(Epetra_DataAccess CV_in, int* Values_in, int Length_in)
  : Epetra_IntSerialDenseMatrix(CV_in, Values_in, Length_in, Length_in, 1)
{
	SetLabel("Epetra::IntSerialDenseVector");
}

//=============================================================================
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(const Epetra_IntSerialDenseVector& Source)
  : Epetra_IntSerialDenseMatrix(Source)
{}

//=============================================================================
Epetra_IntSerialDenseVector::~Epetra_IntSerialDenseVector()
{}

//=========================================================================
Epetra_IntSerialDenseVector& Epetra_IntSerialDenseVector::operator = (const Epetra_IntSerialDenseVector& Source) {
	Epetra_IntSerialDenseMatrix::operator=(Source); // call this->Epetra_IntSerialDenseMatrix::operator =
	return(*this);
}

//=============================================================================
int Epetra_IntSerialDenseVector::MakeViewOf(const Epetra_IntSerialDenseVector& Source) {
	int errorcode = Epetra_IntSerialDenseMatrix::MakeViewOf(Source);
	return(errorcode);
}

//=========================================================================
void Epetra_IntSerialDenseVector::Print(std::ostream& os) const {
	if(CV_ == Copy)
		os << "Data access mode: Copy" << std::endl;
	else
		os << "Data access mode: View" << std::endl;
	if(A_Copied_)
		os << "A_Copied: yes" << std::endl;
	else
		os << "A_Copied: no" << std::endl;
	os << "Length(M): " << M_ << std::endl;
	if(M_ == 0)
		os << "(vector is empty, no values to display)";
	else
		for(int i = 0; i < M_; i++)
      os << (*this)(i) << " ";
	os << std::endl;
}

//=========================================================================
int Epetra_IntSerialDenseVector::Random() {
	int errorcode = Epetra_IntSerialDenseMatrix::Random();
	return(errorcode);
}
