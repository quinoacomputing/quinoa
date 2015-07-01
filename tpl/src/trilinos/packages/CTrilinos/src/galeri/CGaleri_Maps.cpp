
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


#include "CTrilinos_config.h"


#ifdef HAVE_CTRILINOS_GALERI


#include "CTrilinos_enums.h"
#include "CGaleri_Maps.h"
#include "CGaleri_Maps_Cpp.hpp"
#include "Galeri_Maps.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"


//
// Definitions from CGaleri_Maps.h
//


extern "C" {


CT_Epetra_Map_ID_t Galeri_Maps_CreateMap ( 
  char MapType[], CT_Epetra_Comm_ID_t CommID, 
  CT_Teuchos_ParameterList_ID_t ListID )
{
    const Teuchos::RCP<Epetra_Comm> Comm = CEpetra::getComm(CommID);
    const Teuchos::RCP<Teuchos::ParameterList> List = 
        CTeuchos::getParameterList(ListID);
    return CEpetra::storeMap(Galeri::CreateMap(std::string(MapType), *Comm, 
        *List));
}


} // extern "C"


//
// Definitions from CGaleri_Maps_Cpp.hpp
//




#endif /* HAVE_CTRILINOS_GALERI */


