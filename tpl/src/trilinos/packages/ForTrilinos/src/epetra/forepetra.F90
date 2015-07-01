!*********************************************************************
! ForTrilinos: Object-Oriented Fortran 2003 interface to Trilinos
!                Copyright 2010 Sandia Corporation
!
! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
! the U.S. Government retains certain rights in this software.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the Corporation nor the names of the
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov) 
!*********************************************************************


#include "ForTrilinos_config.h"
module forepetra
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  use ForTrilinos_enum_wrappers
  implicit none   ! Prevent implicit typing
#ifdef HAVE_MPI
#include <mpif.h>
#endif

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/epetra/CEpetra*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

!> @name Epetra_Distributor interface
!! @{

  ! _________________ Epetra_Distributor interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Distributor_ID_t Epetra_Distributor_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Distributor_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Distributor_Degeneralize')
    import :: FT_Epetra_Distributor_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Distributor_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Distributor_Generalize ( CT_Epetra_Distributor_ID_t id );

  function Epetra_Distributor_Generalize ( id ) result(that) &
        bind(C,name='Epetra_Distributor_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Distributor_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual Epetra_Distributor * Clone() = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( CT_Epetra_Distributor_ID_t selfID );

  function Epetra_Distributor_Clone ( selfID ) result(that) &
        bind(C,name='Epetra_Distributor_Clone')
    import :: FT_Epetra_Distributor_ID_t
    
    type(FT_Epetra_Distributor_ID_t)                                  :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Distributor();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Distributor_Destroy ( CT_Epetra_Distributor_ID_t * selfID );

  subroutine Epetra_Distributor_Destroy ( selfID ) &
        bind(C,name='Epetra_Distributor_Destroy')
    import :: FT_Epetra_Distributor_ID_t
    
    type(FT_Epetra_Distributor_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual int CreateFromSends( const int & NumExportIDs, const int * ExportPIDs, 
  !!     bool Deterministic, int & NumRemoteIDs ) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_CreateFromSends ( CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  !!     const int * ExportPIDs, boolean Deterministic, int * NumRemoteIDs );

  function Epetra_Distributor_CreateFromSends ( selfID, NumExportIDs, ExportPIDs, &
        Deterministic, NumRemoteIDs ) result(that) &
        bind(C,name='Epetra_Distributor_CreateFromSends')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,FT_boolean_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: NumExportIDs
    integer(c_int)                  ,intent(in)         ,dimension(*) :: ExportPIDs
    integer(FT_boolean_t)           ,intent(in)   ,value              :: Deterministic
    integer(c_int)                  ,intent(inout)                    :: NumRemoteIDs
  end function


  !> <BR> Original C++ prototype:
  !! virtual int CreateFromRecvs( const int & NumRemoteIDs, const int * RemoteGIDs, 
  !!     const int * RemotePIDs, bool Deterministic, int & NumExportIDs, int *& ExportGIDs, 
  !!     int *& ExportPIDs) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_CreateFromRecvs ( CT_Epetra_Distributor_ID_t selfID, int NumRemoteIDs, 
  !!     const int * RemoteGIDs, const int * RemotePIDs, boolean Deterministic, 
  !!     int * NumExportIDs, int ** ExportGIDs, int ** ExportPIDs );

  function Epetra_Distributor_CreateFromRecvs ( selfID, NumRemoteIDs, RemoteGIDs, &
        RemotePIDs, Deterministic, NumExportIDs, ExportGIDs, ExportPIDs ) result(that) &
        bind(C,name='Epetra_Distributor_CreateFromRecvs')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,FT_boolean_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: NumRemoteIDs
    integer(c_int)                  ,intent(in)         ,dimension(*) :: RemoteGIDs
    integer(c_int)                  ,intent(in)         ,dimension(*) :: RemotePIDs
    integer(FT_boolean_t)           ,intent(in)   ,value              :: Deterministic
    integer(c_int)                  ,intent(inout)                    :: NumExportIDs
    integer(c_int)                  ,intent(inout)      ,dimension(*) :: ExportGIDs
    integer(c_int)                  ,intent(inout)      ,dimension(*) :: ExportPIDs
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Do( char * export_objs, int obj_size, int & len_import_objs, 
  !!     char *& import_objs) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_Do ( CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  !!     int obj_size, int * len_import_objs, char ** import_objs );

  function Epetra_Distributor_Do ( selfID, export_objs, obj_size, len_import_objs, &
        import_objs ) result(that) bind(C,name='Epetra_Distributor_Do')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                              ,dimension(*) :: export_objs
    integer(c_int)                  ,intent(in)   ,value              :: obj_size
    integer(c_int)                  ,intent(inout)                    :: len_import_objs
    character(kind=c_char)          ,intent(inout)      ,dimension(*) :: import_objs
  end function


  !> <BR> Original C++ prototype:
  !! virtual int DoReverse( char * export_objs, int obj_size, int & len_import_objs, 
  !!     char *& import_objs ) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_DoReverse ( CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  !!     int obj_size, int * len_import_objs, char ** import_objs );

  function Epetra_Distributor_DoReverse ( selfID, export_objs, obj_size, len_import_objs, &
        import_objs ) result(that) bind(C,name='Epetra_Distributor_DoReverse')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                              ,dimension(*) :: export_objs
    integer(c_int)                  ,intent(in)   ,value              :: obj_size
    integer(c_int)                  ,intent(inout)                    :: len_import_objs
    character(kind=c_char)          ,intent(inout)      ,dimension(*) :: import_objs
  end function


  !> <BR> Original C++ prototype:
  !! virtual int DoPosts( char * export_objs, int obj_size, int & len_import_objs, 
  !!     char *& import_objs ) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_DoPosts ( CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  !!     int obj_size, int * len_import_objs, char ** import_objs );

  function Epetra_Distributor_DoPosts ( selfID, export_objs, obj_size, len_import_objs, &
        import_objs ) result(that) bind(C,name='Epetra_Distributor_DoPosts')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                              ,dimension(*) :: export_objs
    integer(c_int)                  ,intent(in)   ,value              :: obj_size
    integer(c_int)                  ,intent(inout)                    :: len_import_objs
    character(kind=c_char)          ,intent(inout)      ,dimension(*) :: import_objs
  end function


  !> <BR> Original C++ prototype:
  !! virtual int DoWaits() = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID );

  function Epetra_Distributor_DoWaits ( selfID ) result(that) &
        bind(C,name='Epetra_Distributor_DoWaits')
    import :: c_int ,FT_Epetra_Distributor_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int DoReversePosts( char * export_objs, int obj_size, int & len_import_objs, 
  !!     char *& import_objs) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_DoReversePosts ( CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  !!     int obj_size, int * len_import_objs, char ** import_objs );

  function Epetra_Distributor_DoReversePosts ( selfID, export_objs, obj_size, &
        len_import_objs, import_objs ) result(that) &
        bind(C,name='Epetra_Distributor_DoReversePosts')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                              ,dimension(*) :: export_objs
    integer(c_int)                  ,intent(in)   ,value              :: obj_size
    integer(c_int)                  ,intent(inout)                    :: len_import_objs
    character(kind=c_char)          ,intent(inout)      ,dimension(*) :: import_objs
  end function


  !> <BR> Original C++ prototype:
  !! virtual int DoReverseWaits() = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_DoReverseWaits ( CT_Epetra_Distributor_ID_t selfID );

  function Epetra_Distributor_DoReverseWaits ( selfID ) result(that) &
        bind(C,name='Epetra_Distributor_DoReverseWaits')
    import :: c_int ,FT_Epetra_Distributor_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Do( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, 
  !!     char *& import_objs) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_Do_VarLen ( CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  !!     int obj_size, int ** sizes, int * len_import_objs, char ** import_objs );

  function Epetra_Distributor_Do_VarLen ( selfID, export_objs, obj_size, sizes, &
        len_import_objs, import_objs ) result(that) &
        bind(C,name='Epetra_Distributor_Do_VarLen')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                              ,dimension(*) :: export_objs
    integer(c_int)                  ,intent(in)   ,value              :: obj_size
    integer(c_int)                  ,intent(inout)      ,dimension(*) :: sizes
    integer(c_int)                  ,intent(inout)                    :: len_import_objs
    character(kind=c_char)          ,intent(inout)      ,dimension(*) :: import_objs
  end function


  !> <BR> Original C++ prototype:
  !! virtual int DoReverse( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, 
  !!     char *& import_objs) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_DoReverse_VarLen ( CT_Epetra_Distributor_ID_t selfID, 
  !!     char * export_objs, int obj_size, int ** sizes, int * len_import_objs, 
  !!     char ** import_objs );

  function Epetra_Distributor_DoReverse_VarLen ( selfID, export_objs, obj_size, sizes, &
        len_import_objs, import_objs ) result(that) &
        bind(C,name='Epetra_Distributor_DoReverse_VarLen')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                              ,dimension(*) :: export_objs
    integer(c_int)                  ,intent(in)   ,value              :: obj_size
    integer(c_int)                  ,intent(inout)      ,dimension(*) :: sizes
    integer(c_int)                  ,intent(inout)                    :: len_import_objs
    character(kind=c_char)          ,intent(inout)      ,dimension(*) :: import_objs
  end function


  !> <BR> Original C++ prototype:
  !! virtual int DoPosts( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, 
  !!     char *& import_objs) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_DoPosts_VarLen ( CT_Epetra_Distributor_ID_t selfID, char * export_objs, 
  !!     int obj_size, int ** sizes, int * len_import_objs, char ** import_objs );

  function Epetra_Distributor_DoPosts_VarLen ( selfID, export_objs, obj_size, sizes, &
        len_import_objs, import_objs ) result(that) &
        bind(C,name='Epetra_Distributor_DoPosts_VarLen')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                              ,dimension(*) :: export_objs
    integer(c_int)                  ,intent(in)   ,value              :: obj_size
    integer(c_int)                  ,intent(inout)      ,dimension(*) :: sizes
    integer(c_int)                  ,intent(inout)                    :: len_import_objs
    character(kind=c_char)          ,intent(inout)      ,dimension(*) :: import_objs
  end function


  !> <BR> Original C++ prototype:
  !! virtual int DoReversePosts( char * export_objs, int obj_size, int *& sizes, 
  !!     int & len_import_objs, char *& import_objs) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Distributor_DoReversePosts_VarLen ( CT_Epetra_Distributor_ID_t selfID, 
  !!     char * export_objs, int obj_size, int ** sizes, int * len_import_objs, 
  !!     char ** import_objs );

  function Epetra_Distributor_DoReversePosts_VarLen ( selfID, export_objs, obj_size, sizes, &
        len_import_objs, import_objs ) result(that) &
        bind(C,name='Epetra_Distributor_DoReversePosts_VarLen')
    import :: c_int ,FT_Epetra_Distributor_ID_t ,c_char
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_Distributor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                              ,dimension(*) :: export_objs
    integer(c_int)                  ,intent(in)   ,value              :: obj_size
    integer(c_int)                  ,intent(inout)      ,dimension(*) :: sizes
    integer(c_int)                  ,intent(inout)                    :: len_import_objs
    character(kind=c_char)          ,intent(inout)      ,dimension(*) :: import_objs
  end function


!> @}


!> @name Epetra_SerialComm interface
!! @{

  ! _________________ Epetra_SerialComm interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_SerialComm_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_SerialComm_Degeneralize')
    import :: FT_Epetra_SerialComm_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_SerialComm_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_SerialComm_Generalize ( CT_Epetra_SerialComm_ID_t id );

  function Epetra_SerialComm_Generalize ( id ) result(that) &
        bind(C,name='Epetra_SerialComm_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_SerialComm_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialComm();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  );

  function Epetra_SerialComm_Create (  ) result(that) &
        bind(C,name='Epetra_SerialComm_Create')
    import :: FT_Epetra_SerialComm_ID_t
    
    type(FT_Epetra_SerialComm_ID_t)                                  :: that
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialComm(const Epetra_SerialComm& Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( CT_Epetra_SerialComm_ID_t CommID );

  function Epetra_SerialComm_Duplicate ( CommID ) result(that) &
        bind(C,name='Epetra_SerialComm_Duplicate')
    import :: FT_Epetra_SerialComm_ID_t
    
    type(FT_Epetra_SerialComm_ID_t)                                  :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Comm * Clone() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Epetra_SerialComm_Clone ( CT_Epetra_SerialComm_ID_t selfID );

  function Epetra_SerialComm_Clone ( selfID ) result(that) &
        bind(C,name='Epetra_SerialComm_Clone')
    import :: FT_Epetra_Comm_ID_t ,FT_Epetra_SerialComm_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                        :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_SerialComm();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialComm_Destroy ( CT_Epetra_SerialComm_ID_t * selfID );

  subroutine Epetra_SerialComm_Destroy ( selfID ) bind(C,name='Epetra_SerialComm_Destroy')
    import :: FT_Epetra_SerialComm_ID_t
    
    type(FT_Epetra_SerialComm_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void Barrier() const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialComm_Barrier ( CT_Epetra_SerialComm_ID_t selfID );

  subroutine Epetra_SerialComm_Barrier ( selfID ) bind(C,name='Epetra_SerialComm_Barrier')
    import :: FT_Epetra_SerialComm_ID_t
    
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int Broadcast(double * MyVals, int Count, int Root) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_Broadcast_Double ( CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  !!     int Count, int Root );

  function Epetra_SerialComm_Broadcast_Double ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_SerialComm_Broadcast_Double')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_double
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                     ,dimension(*) :: MyVals
    integer(c_int)                 ,intent(in)   ,value              :: Count
    integer(c_int)                 ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! int Broadcast(int * MyVals, int Count, int Root) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_Broadcast_Int ( CT_Epetra_SerialComm_ID_t selfID, int * MyVals, 
  !!     int Count, int Root );

  function Epetra_SerialComm_Broadcast_Int ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_SerialComm_Broadcast_Int')
    import :: c_int ,FT_Epetra_SerialComm_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                     ,dimension(*) :: MyVals
    integer(c_int)                 ,intent(in)   ,value              :: Count
    integer(c_int)                 ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! int Broadcast(long * MyVals, int Count, int Root) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_Broadcast_Long ( CT_Epetra_SerialComm_ID_t selfID, long * MyVals, 
  !!     int Count, int Root );

  function Epetra_SerialComm_Broadcast_Long ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_SerialComm_Broadcast_Long')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_long
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                    ,dimension(*) :: MyVals
    integer(c_int)                 ,intent(in)   ,value              :: Count
    integer(c_int)                 ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! int Broadcast(char * MyVals, int Count, int Root) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_Broadcast_Char ( CT_Epetra_SerialComm_ID_t selfID, char * MyVals, 
  !!     int Count, int Root );

  function Epetra_SerialComm_Broadcast_Char ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_SerialComm_Broadcast_Char')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_char
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                             ,dimension(*) :: MyVals
    integer(c_int)                 ,intent(in)   ,value              :: Count
    integer(c_int)                 ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! int GatherAll(double * MyVals, double * AllVals, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_GatherAll_Double ( CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  !!     double * AllVals, int Count );

  function Epetra_SerialComm_GatherAll_Double ( selfID, MyVals, AllVals, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_GatherAll_Double')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_double
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                     ,dimension(*) :: MyVals
    real(c_double)                                     ,dimension(*) :: AllVals
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int GatherAll(int * MyVals, int * AllVals, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_GatherAll_Int ( CT_Epetra_SerialComm_ID_t selfID, int * MyVals, 
  !!     int * AllVals, int Count );

  function Epetra_SerialComm_GatherAll_Int ( selfID, MyVals, AllVals, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_GatherAll_Int')
    import :: c_int ,FT_Epetra_SerialComm_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                     ,dimension(*) :: MyVals
    integer(c_int)                                     ,dimension(*) :: AllVals
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int GatherAll(long * MyVals, long * AllVals, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_GatherAll_Long ( CT_Epetra_SerialComm_ID_t selfID, long * MyVals, 
  !!     long * AllVals, int Count );

  function Epetra_SerialComm_GatherAll_Long ( selfID, MyVals, AllVals, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_GatherAll_Long')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_long
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                    ,dimension(*) :: MyVals
    integer(c_long)                                    ,dimension(*) :: AllVals
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int SumAll(double * PartialSums, double * GlobalSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_SumAll_Double ( CT_Epetra_SerialComm_ID_t selfID, double * PartialSums, 
  !!     double * GlobalSums, int Count );

  function Epetra_SerialComm_SumAll_Double ( selfID, PartialSums, GlobalSums, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_SumAll_Double')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_double
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                     ,dimension(*) :: PartialSums
    real(c_double)                                     ,dimension(*) :: GlobalSums
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int SumAll(int * PartialSums, int * GlobalSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_SumAll_Int ( CT_Epetra_SerialComm_ID_t selfID, int * PartialSums, 
  !!     int * GlobalSums, int Count );

  function Epetra_SerialComm_SumAll_Int ( selfID, PartialSums, GlobalSums, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_SumAll_Int')
    import :: c_int ,FT_Epetra_SerialComm_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                     ,dimension(*) :: PartialSums
    integer(c_int)                                     ,dimension(*) :: GlobalSums
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int SumAll(long * PartialSums, long * GlobalSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_SumAll_Long ( CT_Epetra_SerialComm_ID_t selfID, long * PartialSums, 
  !!     long * GlobalSums, int Count );

  function Epetra_SerialComm_SumAll_Long ( selfID, PartialSums, GlobalSums, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_SumAll_Long')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_long
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                    ,dimension(*) :: PartialSums
    integer(c_long)                                    ,dimension(*) :: GlobalSums
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_MaxAll_Double ( CT_Epetra_SerialComm_ID_t selfID, double * PartialMaxs, 
  !!     double * GlobalMaxs, int Count );

  function Epetra_SerialComm_MaxAll_Double ( selfID, PartialMaxs, GlobalMaxs, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_MaxAll_Double')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_double
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                     ,dimension(*) :: PartialMaxs
    real(c_double)                                     ,dimension(*) :: GlobalMaxs
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_MaxAll_Int ( CT_Epetra_SerialComm_ID_t selfID, int * PartialMaxs, 
  !!     int * GlobalMaxs, int Count );

  function Epetra_SerialComm_MaxAll_Int ( selfID, PartialMaxs, GlobalMaxs, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_MaxAll_Int')
    import :: c_int ,FT_Epetra_SerialComm_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                     ,dimension(*) :: PartialMaxs
    integer(c_int)                                     ,dimension(*) :: GlobalMaxs
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_MaxAll_Long ( CT_Epetra_SerialComm_ID_t selfID, long * PartialMaxs, 
  !!     long * GlobalMaxs, int Count );

  function Epetra_SerialComm_MaxAll_Long ( selfID, PartialMaxs, GlobalMaxs, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_MaxAll_Long')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_long
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                    ,dimension(*) :: PartialMaxs
    integer(c_long)                                    ,dimension(*) :: GlobalMaxs
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MinAll(double * PartialMins, double * GlobalMins, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_MinAll_Double ( CT_Epetra_SerialComm_ID_t selfID, double * PartialMins, 
  !!     double * GlobalMins, int Count );

  function Epetra_SerialComm_MinAll_Double ( selfID, PartialMins, GlobalMins, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_MinAll_Double')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_double
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                     ,dimension(*) :: PartialMins
    real(c_double)                                     ,dimension(*) :: GlobalMins
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MinAll(int * PartialMins, int * GlobalMins, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_MinAll_Int ( CT_Epetra_SerialComm_ID_t selfID, int * PartialMins, 
  !!     int * GlobalMins, int Count );

  function Epetra_SerialComm_MinAll_Int ( selfID, PartialMins, GlobalMins, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_MinAll_Int')
    import :: c_int ,FT_Epetra_SerialComm_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                     ,dimension(*) :: PartialMins
    integer(c_int)                                     ,dimension(*) :: GlobalMins
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MinAll(long * PartialMins, long * GlobalMins, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_MinAll_Long ( CT_Epetra_SerialComm_ID_t selfID, long * PartialMins, 
  !!     long * GlobalMins, int Count );

  function Epetra_SerialComm_MinAll_Long ( selfID, PartialMins, GlobalMins, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_MinAll_Long')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_long
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                    ,dimension(*) :: PartialMins
    integer(c_long)                                    ,dimension(*) :: GlobalMins
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int ScanSum(double * MyVals, double * ScanSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_ScanSum_Double ( CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  !!     double * ScanSums, int Count );

  function Epetra_SerialComm_ScanSum_Double ( selfID, MyVals, ScanSums, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_ScanSum_Double')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_double
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                     ,dimension(*) :: MyVals
    real(c_double)                                     ,dimension(*) :: ScanSums
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int ScanSum(int * MyVals, int * ScanSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_ScanSum_Int ( CT_Epetra_SerialComm_ID_t selfID, int * MyVals, 
  !!     int * ScanSums, int Count );

  function Epetra_SerialComm_ScanSum_Int ( selfID, MyVals, ScanSums, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_ScanSum_Int')
    import :: c_int ,FT_Epetra_SerialComm_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                     ,dimension(*) :: MyVals
    integer(c_int)                                     ,dimension(*) :: ScanSums
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int ScanSum(long * MyVals, long * ScanSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_ScanSum_Long ( CT_Epetra_SerialComm_ID_t selfID, long * MyVals, 
  !!     long * ScanSums, int Count );

  function Epetra_SerialComm_ScanSum_Long ( selfID, MyVals, ScanSums, Count ) result(that) &
        bind(C,name='Epetra_SerialComm_ScanSum_Long')
    import :: c_int ,FT_Epetra_SerialComm_ID_t ,c_long
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                    ,dimension(*) :: MyVals
    integer(c_long)                                    ,dimension(*) :: ScanSums
    integer(c_int)                 ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MyPID() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_MyPID ( CT_Epetra_SerialComm_ID_t selfID );

  function Epetra_SerialComm_MyPID ( selfID ) result(that) &
        bind(C,name='Epetra_SerialComm_MyPID')
    import :: c_int ,FT_Epetra_SerialComm_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumProc() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialComm_NumProc ( CT_Epetra_SerialComm_ID_t selfID );

  function Epetra_SerialComm_NumProc ( selfID ) result(that) &
        bind(C,name='Epetra_SerialComm_NumProc')
    import :: c_int ,FT_Epetra_SerialComm_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Distributor * CreateDistributor() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Distributor_ID_t Epetra_SerialComm_CreateDistributor ( CT_Epetra_SerialComm_ID_t selfID );

  function Epetra_SerialComm_CreateDistributor ( selfID ) result(that) &
        bind(C,name='Epetra_SerialComm_CreateDistributor')
    import :: FT_Epetra_Distributor_ID_t ,FT_Epetra_SerialComm_ID_t
    
    type(FT_Epetra_Distributor_ID_t)                                  :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Directory_ID_t Epetra_SerialComm_CreateDirectory ( CT_Epetra_SerialComm_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t MapID );

  function Epetra_SerialComm_CreateDirectory ( selfID, MapID ) result(that) &
        bind(C,name='Epetra_SerialComm_CreateDirectory')
    import :: FT_Epetra_Directory_ID_t ,FT_Epetra_SerialComm_ID_t ,FT_Epetra_BlockMap_ID_t
    
    type(FT_Epetra_Directory_ID_t)                                   :: that
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t)  ,intent(in)   ,value              :: MapID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialComm & operator=(const Epetra_SerialComm & Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialComm_Assign ( CT_Epetra_SerialComm_ID_t selfID, 
  !!     CT_Epetra_SerialComm_ID_t CommID );

  subroutine Epetra_SerialComm_Assign ( selfID, CommID ) &
        bind(C,name='Epetra_SerialComm_Assign')
    import :: FT_Epetra_SerialComm_ID_t
    
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialComm_ID_t),intent(in)   ,value              :: CommID
  end subroutine


!> @}


!> @name Epetra_BLAS interface
!! @{

  ! _________________ Epetra_BLAS interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_BLAS_ID_t Epetra_BLAS_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_BLAS_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_BLAS_Degeneralize')
    import :: FT_Epetra_BLAS_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_BLAS_ID_t)                                     :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_BLAS_Generalize ( CT_Epetra_BLAS_ID_t id );

  function Epetra_BLAS_Generalize ( id ) result(that) bind(C,name='Epetra_BLAS_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_BLAS_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_BLAS(void);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BLAS_ID_t Epetra_BLAS_Create (  );

  function Epetra_BLAS_Create (  ) result(that) bind(C,name='Epetra_BLAS_Create')
    import :: FT_Epetra_BLAS_ID_t
    
    type(FT_Epetra_BLAS_ID_t)                                     :: that
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_BLAS(const Epetra_BLAS& BLAS);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BLAS_ID_t Epetra_BLAS_Duplicate ( CT_Epetra_BLAS_ID_t BLASID );

  function Epetra_BLAS_Duplicate ( BLASID ) result(that) &
        bind(C,name='Epetra_BLAS_Duplicate')
    import :: FT_Epetra_BLAS_ID_t
    
    type(FT_Epetra_BLAS_ID_t)                                     :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: BLASID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_BLAS(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_Destroy ( CT_Epetra_BLAS_ID_t * selfID );

  subroutine Epetra_BLAS_Destroy ( selfID ) bind(C,name='Epetra_BLAS_Destroy')
    import :: FT_Epetra_BLAS_ID_t
    
    type(FT_Epetra_BLAS_ID_t)                                     :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! float ASUM(const int N, const float * X, const int INCX = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! float Epetra_BLAS_ASUM_Float ( CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  !!     const int INCX );

  function Epetra_BLAS_ASUM_Float ( selfID, N, X, INCX ) result(that) &
        bind(C,name='Epetra_BLAS_ASUM_Float')
    import :: c_float ,FT_Epetra_BLAS_ID_t ,c_int
    
    real(c_float)                                                 :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: INCX
  end function


  !> <BR> Original C++ prototype:
  !! double ASUM(const int N, const double * X, const int INCX = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_BLAS_ASUM_Double ( CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  !!     const int INCX );

  function Epetra_BLAS_ASUM_Double ( selfID, N, X, INCX ) result(that) &
        bind(C,name='Epetra_BLAS_ASUM_Double')
    import :: c_double ,FT_Epetra_BLAS_ID_t ,c_int
    
    real(c_double)                                                :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: INCX
  end function


  !> <BR> Original C++ prototype:
  !! float DOT(const int N, const float * X, const float * Y, const int INCX = 1, const int INCY = 
  !!     1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! float Epetra_BLAS_DOT_Float ( CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  !!     const float * Y, const int INCX, const int INCY );

  function Epetra_BLAS_DOT_Float ( selfID, N, X, Y, INCX, INCY ) result(that) &
        bind(C,name='Epetra_BLAS_DOT_Float')
    import :: c_float ,FT_Epetra_BLAS_ID_t ,c_int
    
    real(c_float)                                                 :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: X
    real(c_float)               ,intent(in)         ,dimension(*) :: Y
    integer(c_int)              ,intent(in)   ,value              :: INCX
    integer(c_int)              ,intent(in)   ,value              :: INCY
  end function


  !> <BR> Original C++ prototype:
  !! double DOT(const int N, const double * X, const double * Y, const int INCX = 1, 
  !!     const int INCY = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_BLAS_DOT_Double ( CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  !!     const double * Y, const int INCX, const int INCY );

  function Epetra_BLAS_DOT_Double ( selfID, N, X, Y, INCX, INCY ) result(that) &
        bind(C,name='Epetra_BLAS_DOT_Double')
    import :: c_double ,FT_Epetra_BLAS_ID_t ,c_int
    
    real(c_double)                                                :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: X
    real(c_double)              ,intent(in)         ,dimension(*) :: Y
    integer(c_int)              ,intent(in)   ,value              :: INCX
    integer(c_int)              ,intent(in)   ,value              :: INCY
  end function


  !> <BR> Original C++ prototype:
  !! float NRM2(const int N, const float * X, const int INCX = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! float Epetra_BLAS_NRM2_Float ( CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  !!     const int INCX );

  function Epetra_BLAS_NRM2_Float ( selfID, N, X, INCX ) result(that) &
        bind(C,name='Epetra_BLAS_NRM2_Float')
    import :: c_float ,FT_Epetra_BLAS_ID_t ,c_int
    
    real(c_float)                                                 :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: INCX
  end function


  !> <BR> Original C++ prototype:
  !! double NRM2(const int N, const double * X, const int INCX = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_BLAS_NRM2_Double ( CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  !!     const int INCX );

  function Epetra_BLAS_NRM2_Double ( selfID, N, X, INCX ) result(that) &
        bind(C,name='Epetra_BLAS_NRM2_Double')
    import :: c_double ,FT_Epetra_BLAS_ID_t ,c_int
    
    real(c_double)                                                :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: INCX
  end function


  !> <BR> Original C++ prototype:
  !! void SCAL( const int N, const float ALPHA, float * X, const int INCX = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_SCAL_Float ( CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  !!     float * X, const int INCX );

  subroutine Epetra_BLAS_SCAL_Float ( selfID, N, ALPHA, X, INCX ) &
        bind(C,name='Epetra_BLAS_SCAL_Float')
    import :: FT_Epetra_BLAS_ID_t ,c_int ,c_float
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)   ,value              :: ALPHA
    real(c_float)                                   ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: INCX
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SCAL( const int N, const double ALPHA, double * X, const int INCX = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_SCAL_Double ( CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  !!     double * X, const int INCX );

  subroutine Epetra_BLAS_SCAL_Double ( selfID, N, ALPHA, X, INCX ) &
        bind(C,name='Epetra_BLAS_SCAL_Double')
    import :: FT_Epetra_BLAS_ID_t ,c_int ,c_double
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)   ,value              :: ALPHA
    real(c_double)                                  ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: INCX
  end subroutine


  !> <BR> Original C++ prototype:
  !! void COPY( const int N, const float * X, float * Y, const int INCX = 1, const int INCY = 
  !!     1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_COPY_Float ( CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  !!     float * Y, const int INCX, const int INCY );

  subroutine Epetra_BLAS_COPY_Float ( selfID, N, X, Y, INCX, INCY ) &
        bind(C,name='Epetra_BLAS_COPY_Float')
    import :: FT_Epetra_BLAS_ID_t ,c_int ,c_float
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: X
    real(c_float)                                   ,dimension(*) :: Y
    integer(c_int)              ,intent(in)   ,value              :: INCX
    integer(c_int)              ,intent(in)   ,value              :: INCY
  end subroutine


  !> <BR> Original C++ prototype:
  !! void COPY( const int N, const double * X, double * Y, const int INCX = 1, const int INCY = 
  !!     1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_COPY_Double ( CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  !!     double * Y, const int INCX, const int INCY );

  subroutine Epetra_BLAS_COPY_Double ( selfID, N, X, Y, INCX, INCY ) &
        bind(C,name='Epetra_BLAS_COPY_Double')
    import :: FT_Epetra_BLAS_ID_t ,c_int ,c_double
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: X
    real(c_double)                                  ,dimension(*) :: Y
    integer(c_int)              ,intent(in)   ,value              :: INCX
    integer(c_int)              ,intent(in)   ,value              :: INCY
  end subroutine


  !> <BR> Original C++ prototype:
  !! int IAMAX( const int N, const float * X, const int INCX = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BLAS_IAMAX_Float ( CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  !!     const int INCX );

  function Epetra_BLAS_IAMAX_Float ( selfID, N, X, INCX ) result(that) &
        bind(C,name='Epetra_BLAS_IAMAX_Float')
    import :: c_int ,FT_Epetra_BLAS_ID_t ,c_float
    
    integer(c_int)                                                :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: INCX
  end function


  !> <BR> Original C++ prototype:
  !! int IAMAX( const int N, const double * X, const int INCX = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BLAS_IAMAX_Double ( CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  !!     const int INCX );

  function Epetra_BLAS_IAMAX_Double ( selfID, N, X, INCX ) result(that) &
        bind(C,name='Epetra_BLAS_IAMAX_Double')
    import :: c_int ,FT_Epetra_BLAS_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: INCX
  end function


  !> <BR> Original C++ prototype:
  !! void AXPY( const int N, const float ALPHA, const float * X, float * Y, const int INCX = 1, 
  !!     const int INCY = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_AXPY_Float ( CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  !!     const float * X, float * Y, const int INCX, const int INCY );

  subroutine Epetra_BLAS_AXPY_Float ( selfID, N, ALPHA, X, Y, INCX, INCY ) &
        bind(C,name='Epetra_BLAS_AXPY_Float')
    import :: FT_Epetra_BLAS_ID_t ,c_int ,c_float
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)   ,value              :: ALPHA
    real(c_float)               ,intent(in)         ,dimension(*) :: X
    real(c_float)                                   ,dimension(*) :: Y
    integer(c_int)              ,intent(in)   ,value              :: INCX
    integer(c_int)              ,intent(in)   ,value              :: INCY
  end subroutine


  !> <BR> Original C++ prototype:
  !! void AXPY( const int N, const double ALPHA, const double * X, double * Y, const int INCX = 1, 
  !!     const int INCY = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_AXPY_Double ( CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  !!     const double * X, double * Y, const int INCX, const int INCY );

  subroutine Epetra_BLAS_AXPY_Double ( selfID, N, ALPHA, X, Y, INCX, INCY ) &
        bind(C,name='Epetra_BLAS_AXPY_Double')
    import :: FT_Epetra_BLAS_ID_t ,c_int ,c_double
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)   ,value              :: ALPHA
    real(c_double)              ,intent(in)         ,dimension(*) :: X
    real(c_double)                                  ,dimension(*) :: Y
    integer(c_int)              ,intent(in)   ,value              :: INCX
    integer(c_int)              ,intent(in)   ,value              :: INCY
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEMV(const char TRANS, const int M, const int N, const float ALPHA, const float * A, 
  !!     const int LDA, const float * X, const float BETA, float * Y, const int INCX = 1, 
  !!     const int INCY = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_GEMV_Float ( CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  !!     const int N, const float ALPHA, const float * A, const int LDA, const float * X, 
  !!     const float BETA, float * Y, const int INCX, const int INCY );

  subroutine Epetra_BLAS_GEMV_Float ( selfID, TRANS, M, N, ALPHA, A, LDA, X, BETA, Y, INCX, &
        INCY ) bind(C,name='Epetra_BLAS_GEMV_Float')
    import :: FT_Epetra_BLAS_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)   ,value              :: ALPHA
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)         ,dimension(*) :: X
    real(c_float)               ,intent(in)   ,value              :: BETA
    real(c_float)                                   ,dimension(*) :: Y
    integer(c_int)              ,intent(in)   ,value              :: INCX
    integer(c_int)              ,intent(in)   ,value              :: INCY
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEMV(const char TRANS, const int M, const int N, const double ALPHA, const double * A, 
  !!     const int LDA, const double * X, const double BETA, double * Y, const int INCX = 1, 
  !!     const int INCY = 1) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_GEMV_Double ( CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  !!     const int N, const double ALPHA, const double * A, const int LDA, const double * X, 
  !!     const double BETA, double * Y, const int INCX, const int INCY );

  subroutine Epetra_BLAS_GEMV_Double ( selfID, TRANS, M, N, ALPHA, A, LDA, X, BETA, Y, INCX, &
        INCY ) bind(C,name='Epetra_BLAS_GEMV_Double')
    import :: FT_Epetra_BLAS_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)   ,value              :: ALPHA
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)         ,dimension(*) :: X
    real(c_double)              ,intent(in)   ,value              :: BETA
    real(c_double)                                  ,dimension(*) :: Y
    integer(c_int)              ,intent(in)   ,value              :: INCX
    integer(c_int)              ,intent(in)   ,value              :: INCY
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEMM(const char TRANSA, const char TRANSB, const int M, const int N, const int K, 
  !!     const float ALPHA, const float * A, const int LDA, const float * B, const int LDB, 
  !!     const float BETA, float * C, const int LDC) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_GEMM_Float ( CT_Epetra_BLAS_ID_t selfID, const char TRANSA, 
  !!     const char TRANSB, const int M, const int N, const int K, const float ALPHA, 
  !!     const float * A, const int LDA, const float * B, const int LDB, const float BETA, 
  !!     float * C, const int LDC );

  subroutine Epetra_BLAS_GEMM_Float ( selfID, TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, &
        LDB, BETA, C, LDC ) bind(C,name='Epetra_BLAS_GEMM_Float')
    import :: FT_Epetra_BLAS_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: TRANSA
    character(kind=c_char)      ,intent(in)   ,value              :: TRANSB
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: K
    real(c_float)               ,intent(in)   ,value              :: ALPHA
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)         ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)               ,intent(in)   ,value              :: BETA
    real(c_float)                                   ,dimension(*) :: C
    integer(c_int)              ,intent(in)   ,value              :: LDC
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEMM(const char TRANSA, const char TRANSB, const int M, const int N, const int K, 
  !!     const double ALPHA, const double * A, const int LDA, const double * B, const int LDB, 
  !!     const double BETA, double * C, const int LDC) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_GEMM_Double ( CT_Epetra_BLAS_ID_t selfID, const char TRANSA, 
  !!     const char TRANSB, const int M, const int N, const int K, const double ALPHA, 
  !!     const double * A, const int LDA, const double * B, const int LDB, const double BETA, 
  !!     double * C, const int LDC );

  subroutine Epetra_BLAS_GEMM_Double ( selfID, TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, &
        LDB, BETA, C, LDC ) bind(C,name='Epetra_BLAS_GEMM_Double')
    import :: FT_Epetra_BLAS_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: TRANSA
    character(kind=c_char)      ,intent(in)   ,value              :: TRANSB
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: K
    real(c_double)              ,intent(in)   ,value              :: ALPHA
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)         ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)              ,intent(in)   ,value              :: BETA
    real(c_double)                                  ,dimension(*) :: C
    integer(c_int)              ,intent(in)   ,value              :: LDC
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYMM(const char SIDE, const char UPLO, const int M, const int N, const float ALPHA, 
  !!     const float * A, const int LDA, const float * B, const int LDB, const float BETA, 
  !!     float * C, const int LDC) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_SYMM_Float ( CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  !!     const int M, const int N, const float ALPHA, const float * A, const int LDA, 
  !!     const float * B, const int LDB, const float BETA, float * C, const int LDC );

  subroutine Epetra_BLAS_SYMM_Float ( selfID, SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, &
        C, LDC ) bind(C,name='Epetra_BLAS_SYMM_Float')
    import :: FT_Epetra_BLAS_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: SIDE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)   ,value              :: ALPHA
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)         ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)               ,intent(in)   ,value              :: BETA
    real(c_float)                                   ,dimension(*) :: C
    integer(c_int)              ,intent(in)   ,value              :: LDC
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYMM(const char SIDE, const char UPLO, const int M, const int N, const double ALPHA, 
  !!     const double * A, const int LDA, const double * B, const int LDB, const double BETA, 
  !!     double * C, const int LDC) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_SYMM_Double ( CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  !!     const int M, const int N, const double ALPHA, const double * A, const int LDA, 
  !!     const double * B, const int LDB, const double BETA, double * C, const int LDC );

  subroutine Epetra_BLAS_SYMM_Double ( selfID, SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, &
        BETA, C, LDC ) bind(C,name='Epetra_BLAS_SYMM_Double')
    import :: FT_Epetra_BLAS_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: SIDE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)   ,value              :: ALPHA
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)         ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)              ,intent(in)   ,value              :: BETA
    real(c_double)                                  ,dimension(*) :: C
    integer(c_int)              ,intent(in)   ,value              :: LDC
  end subroutine


  !> <BR> Original C++ prototype:
  !! void TRMM(const char SIDE, const char UPLO, const char TRANSA, const char DIAG, const int M, 
  !!     const int N, const float ALPHA, const float * A, const int LDA, float * B, 
  !!     const int LDB) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_TRMM_Float ( CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  !!     const char TRANSA, const char DIAG, const int M, const int N, const float ALPHA, 
  !!     const float * A, const int LDA, float * B, const int LDB );

  subroutine Epetra_BLAS_TRMM_Float ( selfID, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
        B, LDB ) bind(C,name='Epetra_BLAS_TRMM_Float')
    import :: FT_Epetra_BLAS_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: SIDE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    character(kind=c_char)      ,intent(in)   ,value              :: TRANSA
    character(kind=c_char)      ,intent(in)   ,value              :: DIAG
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)   ,value              :: ALPHA
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
  end subroutine


  !> <BR> Original C++ prototype:
  !! void TRMM(const char SIDE, const char UPLO, const char TRANSA, const char DIAG, const int M, 
  !!     const int N, const double ALPHA, const double * A, const int LDA, double * B, 
  !!     const int LDB) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BLAS_TRMM_Double ( CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  !!     const char TRANSA, const char DIAG, const int M, const int N, const double ALPHA, 
  !!     const double * A, const int LDA, double * B, const int LDB );

  subroutine Epetra_BLAS_TRMM_Double ( selfID, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, &
        LDA, B, LDB ) bind(C,name='Epetra_BLAS_TRMM_Double')
    import :: FT_Epetra_BLAS_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_BLAS_ID_t)   ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: SIDE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    character(kind=c_char)      ,intent(in)   ,value              :: TRANSA
    character(kind=c_char)      ,intent(in)   ,value              :: DIAG
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)   ,value              :: ALPHA
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
  end subroutine


!> @}


!> @name Epetra_Comm interface
!! @{

  ! _________________ Epetra_Comm interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Epetra_Comm_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Comm_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Comm_Degeneralize')
    import :: FT_Epetra_Comm_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                     :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Comm_Generalize ( CT_Epetra_Comm_ID_t id );

  function Epetra_Comm_Generalize ( id ) result(that) bind(C,name='Epetra_Comm_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Comm_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual Epetra_Comm * Clone() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID );

  function Epetra_Comm_Clone ( selfID ) result(that) bind(C,name='Epetra_Comm_Clone')
    import :: FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                     :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Comm();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID );

  subroutine Epetra_Comm_Destroy ( selfID ) bind(C,name='Epetra_Comm_Destroy')
    import :: FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                     :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual void Barrier() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Comm_Barrier ( CT_Epetra_Comm_ID_t selfID );

  subroutine Epetra_Comm_Barrier ( selfID ) bind(C,name='Epetra_Comm_Barrier')
    import :: FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual int Broadcast(double * MyVals, int Count, int Root) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_Broadcast_Double ( CT_Epetra_Comm_ID_t selfID, double * MyVals, int Count, 
  !!     int Root );

  function Epetra_Comm_Broadcast_Double ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_Comm_Broadcast_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: MyVals
    integer(c_int)              ,intent(in)   ,value              :: Count
    integer(c_int)              ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Broadcast(int * MyVals, int Count, int Root) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_Broadcast_Int ( CT_Epetra_Comm_ID_t selfID, int * MyVals, int Count, 
  !!     int Root );

  function Epetra_Comm_Broadcast_Int ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_Comm_Broadcast_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: MyVals
    integer(c_int)              ,intent(in)   ,value              :: Count
    integer(c_int)              ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Broadcast(long * MyVals, int Count, int Root) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_Broadcast_Long ( CT_Epetra_Comm_ID_t selfID, long * MyVals, int Count, 
  !!     int Root );

  function Epetra_Comm_Broadcast_Long ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_Comm_Broadcast_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: MyVals
    integer(c_int)              ,intent(in)   ,value              :: Count
    integer(c_int)              ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Broadcast(char * MyVals, int Count, int Root) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_Broadcast_Char ( CT_Epetra_Comm_ID_t selfID, char * MyVals, int Count, 
  !!     int Root );

  function Epetra_Comm_Broadcast_Char ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_Comm_Broadcast_Char')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_char
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    character(kind=c_char)                          ,dimension(*) :: MyVals
    integer(c_int)              ,intent(in)   ,value              :: Count
    integer(c_int)              ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! virtual int GatherAll(double * MyVals, double * AllVals, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_GatherAll_Double ( CT_Epetra_Comm_ID_t selfID, double * MyVals, 
  !!     double * AllVals, int Count );

  function Epetra_Comm_GatherAll_Double ( selfID, MyVals, AllVals, Count ) result(that) &
        bind(C,name='Epetra_Comm_GatherAll_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: MyVals
    real(c_double)                                  ,dimension(*) :: AllVals
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int GatherAll(int * MyVals, int * AllVals, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_GatherAll_Int ( CT_Epetra_Comm_ID_t selfID, int * MyVals, int * AllVals, 
  !!     int Count );

  function Epetra_Comm_GatherAll_Int ( selfID, MyVals, AllVals, Count ) result(that) &
        bind(C,name='Epetra_Comm_GatherAll_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: MyVals
    integer(c_int)                                  ,dimension(*) :: AllVals
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int GatherAll(long * MyVals, long * AllVals, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_GatherAll_Long ( CT_Epetra_Comm_ID_t selfID, long * MyVals, long * AllVals, 
  !!     int Count );

  function Epetra_Comm_GatherAll_Long ( selfID, MyVals, AllVals, Count ) result(that) &
        bind(C,name='Epetra_Comm_GatherAll_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: MyVals
    integer(c_long)                                 ,dimension(*) :: AllVals
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int SumAll(double * PartialSums, double * GlobalSums, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_SumAll_Double ( CT_Epetra_Comm_ID_t selfID, double * PartialSums, 
  !!     double * GlobalSums, int Count );

  function Epetra_Comm_SumAll_Double ( selfID, PartialSums, GlobalSums, Count ) result(that) &
        bind(C,name='Epetra_Comm_SumAll_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: PartialSums
    real(c_double)                                  ,dimension(*) :: GlobalSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int SumAll(int * PartialSums, int * GlobalSums, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_SumAll_Int ( CT_Epetra_Comm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  !!     int Count );

  function Epetra_Comm_SumAll_Int ( selfID, PartialSums, GlobalSums, Count ) result(that) &
        bind(C,name='Epetra_Comm_SumAll_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: PartialSums
    integer(c_int)                                  ,dimension(*) :: GlobalSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int SumAll(long * PartialSums, long * GlobalSums, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_SumAll_Long ( CT_Epetra_Comm_ID_t selfID, long * PartialSums, 
  !!     long * GlobalSums, int Count );

  function Epetra_Comm_SumAll_Long ( selfID, PartialSums, GlobalSums, Count ) result(that) &
        bind(C,name='Epetra_Comm_SumAll_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: PartialSums
    integer(c_long)                                 ,dimension(*) :: GlobalSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_MaxAll_Double ( CT_Epetra_Comm_ID_t selfID, double * PartialMaxs, 
  !!     double * GlobalMaxs, int Count );

  function Epetra_Comm_MaxAll_Double ( selfID, PartialMaxs, GlobalMaxs, Count ) result(that) &
        bind(C,name='Epetra_Comm_MaxAll_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: PartialMaxs
    real(c_double)                                  ,dimension(*) :: GlobalMaxs
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_MaxAll_Int ( CT_Epetra_Comm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  !!     int Count );

  function Epetra_Comm_MaxAll_Int ( selfID, PartialMaxs, GlobalMaxs, Count ) result(that) &
        bind(C,name='Epetra_Comm_MaxAll_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: PartialMaxs
    integer(c_int)                                  ,dimension(*) :: GlobalMaxs
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_MaxAll_Long ( CT_Epetra_Comm_ID_t selfID, long * PartialMaxs, 
  !!     long * GlobalMaxs, int Count );

  function Epetra_Comm_MaxAll_Long ( selfID, PartialMaxs, GlobalMaxs, Count ) result(that) &
        bind(C,name='Epetra_Comm_MaxAll_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: PartialMaxs
    integer(c_long)                                 ,dimension(*) :: GlobalMaxs
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int MinAll(double * PartialMins, double * GlobalMins, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_MinAll_Double ( CT_Epetra_Comm_ID_t selfID, double * PartialMins, 
  !!     double * GlobalMins, int Count );

  function Epetra_Comm_MinAll_Double ( selfID, PartialMins, GlobalMins, Count ) result(that) &
        bind(C,name='Epetra_Comm_MinAll_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: PartialMins
    real(c_double)                                  ,dimension(*) :: GlobalMins
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int MinAll(int * PartialMins, int * GlobalMins, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_MinAll_Int ( CT_Epetra_Comm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  !!     int Count );

  function Epetra_Comm_MinAll_Int ( selfID, PartialMins, GlobalMins, Count ) result(that) &
        bind(C,name='Epetra_Comm_MinAll_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: PartialMins
    integer(c_int)                                  ,dimension(*) :: GlobalMins
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int MinAll(long * PartialMins, long * GlobalMins, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_MinAll_Long ( CT_Epetra_Comm_ID_t selfID, long * PartialMins, 
  !!     long * GlobalMins, int Count );

  function Epetra_Comm_MinAll_Long ( selfID, PartialMins, GlobalMins, Count ) result(that) &
        bind(C,name='Epetra_Comm_MinAll_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: PartialMins
    integer(c_long)                                 ,dimension(*) :: GlobalMins
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ScanSum(double * MyVals, double * ScanSums, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_ScanSum_Double ( CT_Epetra_Comm_ID_t selfID, double * MyVals, 
  !!     double * ScanSums, int Count );

  function Epetra_Comm_ScanSum_Double ( selfID, MyVals, ScanSums, Count ) result(that) &
        bind(C,name='Epetra_Comm_ScanSum_Double')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: MyVals
    real(c_double)                                  ,dimension(*) :: ScanSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ScanSum(int * MyVals, int * ScanSums, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_ScanSum_Int ( CT_Epetra_Comm_ID_t selfID, int * MyVals, int * ScanSums, 
  !!     int Count );

  function Epetra_Comm_ScanSum_Int ( selfID, MyVals, ScanSums, Count ) result(that) &
        bind(C,name='Epetra_Comm_ScanSum_Int')
    import :: c_int ,FT_Epetra_Comm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: MyVals
    integer(c_int)                                  ,dimension(*) :: ScanSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ScanSum(long * MyVals, long * ScanSums, int Count) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_ScanSum_Long ( CT_Epetra_Comm_ID_t selfID, long * MyVals, long * ScanSums, 
  !!     int Count );

  function Epetra_Comm_ScanSum_Long ( selfID, MyVals, ScanSums, Count ) result(that) &
        bind(C,name='Epetra_Comm_ScanSum_Long')
    import :: c_int ,FT_Epetra_Comm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: MyVals
    integer(c_long)                                 ,dimension(*) :: ScanSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! virtual int MyPID() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_MyPID ( CT_Epetra_Comm_ID_t selfID );

  function Epetra_Comm_MyPID ( selfID ) result(that) bind(C,name='Epetra_Comm_MyPID')
    import :: c_int ,FT_Epetra_Comm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumProc() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Comm_NumProc ( CT_Epetra_Comm_ID_t selfID );

  function Epetra_Comm_NumProc ( selfID ) result(that) bind(C,name='Epetra_Comm_NumProc')
    import :: c_int ,FT_Epetra_Comm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual Epetra_Distributor * CreateDistributor() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Distributor_ID_t Epetra_Comm_CreateDistributor ( CT_Epetra_Comm_ID_t selfID );

  function Epetra_Comm_CreateDistributor ( selfID ) result(that) &
        bind(C,name='Epetra_Comm_CreateDistributor')
    import :: FT_Epetra_Distributor_ID_t ,FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_Distributor_ID_t)                                  :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Directory_ID_t Epetra_Comm_CreateDirectory ( CT_Epetra_Comm_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t MapID );

  function Epetra_Comm_CreateDirectory ( selfID, MapID ) result(that) &
        bind(C,name='Epetra_Comm_CreateDirectory')
    import :: FT_Epetra_Directory_ID_t ,FT_Epetra_Comm_ID_t ,FT_Epetra_BlockMap_ID_t
    
    type(FT_Epetra_Directory_ID_t)                                  :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: MapID
  end function


!> @}


!> @name Epetra_Operator interface
!! @{

  ! _________________ Epetra_Operator interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Operator_ID_t Epetra_Operator_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Operator_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Operator_Degeneralize')
    import :: FT_Epetra_Operator_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Operator_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Operator_Generalize ( CT_Epetra_Operator_ID_t id );

  function Epetra_Operator_Generalize ( id ) result(that) &
        bind(C,name='Epetra_Operator_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Operator_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Operator();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Operator_Destroy ( CT_Epetra_Operator_ID_t * selfID );

  subroutine Epetra_Operator_Destroy ( selfID ) bind(C,name='Epetra_Operator_Destroy')
    import :: FT_Epetra_Operator_ID_t
    
    type(FT_Epetra_Operator_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual int SetUseTranspose(bool UseTranspose) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Operator_SetUseTranspose ( CT_Epetra_Operator_ID_t selfID, boolean UseTranspose );

  function Epetra_Operator_SetUseTranspose ( selfID, UseTranspose ) result(that) &
        bind(C,name='Epetra_Operator_SetUseTranspose')
    import :: c_int ,FT_Epetra_Operator_ID_t ,FT_boolean_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)        ,intent(in)   ,value              :: UseTranspose
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Operator_Apply ( CT_Epetra_Operator_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  !!     CT_Epetra_MultiVector_ID_t YID );

  function Epetra_Operator_Apply ( selfID, XID, YID ) result(that) &
        bind(C,name='Epetra_Operator_Apply')
    import :: c_int ,FT_Epetra_Operator_ID_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Operator_ApplyInverse ( CT_Epetra_Operator_ID_t selfID, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

  function Epetra_Operator_ApplyInverse ( selfID, XID, YID ) result(that) &
        bind(C,name='Epetra_Operator_ApplyInverse')
    import :: c_int ,FT_Epetra_Operator_ID_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double NormInf() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_Operator_NormInf ( CT_Epetra_Operator_ID_t selfID );

  function Epetra_Operator_NormInf ( selfID ) result(that) &
        bind(C,name='Epetra_Operator_NormInf')
    import :: c_double ,FT_Epetra_Operator_ID_t
    
    real(c_double)                                                 :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const char * Label() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Epetra_Operator_Label ( CT_Epetra_Operator_ID_t selfID );

  function Epetra_Operator_Label ( selfID ) result(that) &
        bind(C,name='Epetra_Operator_Label')
    import :: c_ptr ,FT_Epetra_Operator_ID_t
    
    type(c_ptr)                                                    :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool UseTranspose() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_Operator_UseTranspose ( CT_Epetra_Operator_ID_t selfID );

  function Epetra_Operator_UseTranspose ( selfID ) result(that) &
        bind(C,name='Epetra_Operator_UseTranspose')
    import :: FT_boolean_t ,FT_Epetra_Operator_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool HasNormInf() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_Operator_HasNormInf ( CT_Epetra_Operator_ID_t selfID );

  function Epetra_Operator_HasNormInf ( selfID ) result(that) &
        bind(C,name='Epetra_Operator_HasNormInf')
    import :: FT_boolean_t ,FT_Epetra_Operator_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_Comm & Comm() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Epetra_Operator_Comm ( CT_Epetra_Operator_ID_t selfID );

  function Epetra_Operator_Comm ( selfID ) result(that) bind(C,name='Epetra_Operator_Comm')
    import :: FT_Epetra_Comm_ID_t ,FT_Epetra_Operator_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                      :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_Map & OperatorDomainMap() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_Operator_OperatorDomainMap ( CT_Epetra_Operator_ID_t selfID );

  function Epetra_Operator_OperatorDomainMap ( selfID ) result(that) &
        bind(C,name='Epetra_Operator_OperatorDomainMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_Operator_ID_t
    
    type(FT_Epetra_Map_ID_t)                                       :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_Map & OperatorRangeMap() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_Operator_OperatorRangeMap ( CT_Epetra_Operator_ID_t selfID );

  function Epetra_Operator_OperatorRangeMap ( selfID ) result(that) &
        bind(C,name='Epetra_Operator_OperatorRangeMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_Operator_ID_t
    
    type(FT_Epetra_Map_ID_t)                                       :: that
    type(FT_Epetra_Operator_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_MultiVector interface
!! @{

  ! _________________ Epetra_MultiVector interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_MultiVector_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_MultiVector_Degeneralize')
    import :: FT_Epetra_MultiVector_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_MultiVector_Generalize ( CT_Epetra_MultiVector_ID_t id );

  function Epetra_MultiVector_Generalize ( id ) result(that) &
        bind(C,name='Epetra_MultiVector_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut = true);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( CT_Epetra_BlockMap_ID_t MapID, 
  !!     int NumVectors, boolean zeroOut );

  function Epetra_MultiVector_Create ( MapID, NumVectors, zeroOut ) result(that) &
        bind(C,name='Epetra_MultiVector_Create')
    import :: FT_Epetra_MultiVector_ID_t ,FT_Epetra_BlockMap_ID_t ,c_int ,FT_boolean_t
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    type(FT_Epetra_BlockMap_ID_t)   ,intent(in)   ,value              :: MapID
    integer(c_int)                  ,intent(in)   ,value              :: NumVectors
    integer(FT_boolean_t)           ,intent(in)   ,value              :: zeroOut
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector(const Epetra_MultiVector& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( CT_Epetra_MultiVector_ID_t SourceID );

  function Epetra_MultiVector_Duplicate ( SourceID ) result(that) &
        bind(C,name='Epetra_MultiVector_Duplicate')
    import :: FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: SourceID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double *A, int MyLDA, 
  !!     int NumVectors);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_From2DA ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_BlockMap_ID_t MapID, double * A, int MyLDA, int NumVectors );

  function Epetra_MultiVector_Create_From2DA ( CV, MapID, A, MyLDA, NumVectors ) result(that) &
        bind(C,name='Epetra_MultiVector_Create_From2DA')
    import :: FT_Epetra_MultiVector_ID_t ,FT_Epetra_DataAccess_E_t , &
          FT_Epetra_BlockMap_ID_t ,c_double ,c_int
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_BlockMap_ID_t)   ,intent(in)   ,value              :: MapID
    real(c_double)                                      ,dimension(*) :: A
    integer(c_int)                  ,intent(in)   ,value              :: MyLDA
    integer(c_int)                  ,intent(in)   ,value              :: NumVectors
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double **ArrayOfPointers, 
  !!     int NumVectors);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromAOP ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_BlockMap_ID_t MapID, double ** ArrayOfPointers, int NumVectors );

  function Epetra_MultiVector_Create_FromAOP ( CV, MapID, ArrayOfPointers, NumVectors ) result(that) &
        bind(C,name='Epetra_MultiVector_Create_FromAOP')
    import :: FT_Epetra_MultiVector_ID_t ,FT_Epetra_DataAccess_E_t , &
          FT_Epetra_BlockMap_ID_t ,c_double ,c_int
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_BlockMap_ID_t)   ,intent(in)   ,value              :: MapID
    real(c_double)                                      ,dimension(*) :: ArrayOfPointers
    integer(c_int)                  ,intent(in)   ,value              :: NumVectors
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int *Indices, 
  !!     int NumVectors);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromList ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_MultiVector_ID_t SourceID, int * Indices, int NumVectors );

  function Epetra_MultiVector_Create_FromList ( CV, SourceID, Indices, NumVectors ) result(that) &
        bind(C,name='Epetra_MultiVector_Create_FromList')
    import :: FT_Epetra_MultiVector_ID_t ,FT_Epetra_DataAccess_E_t ,c_int
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: SourceID
    integer(c_int)                                      ,dimension(*) :: Indices
    integer(c_int)                  ,intent(in)   ,value              :: NumVectors
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int StartIndex, 
  !!     int NumVectors);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromRange ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_MultiVector_ID_t SourceID, int StartIndex, int NumVectors );

  function Epetra_MultiVector_Create_FromRange ( CV, SourceID, StartIndex, NumVectors ) result(that) &
        bind(C,name='Epetra_MultiVector_Create_FromRange')
    import :: FT_Epetra_MultiVector_ID_t ,FT_Epetra_DataAccess_E_t ,c_int
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: SourceID
    integer(c_int)                  ,intent(in)   ,value              :: StartIndex
    integer(c_int)                  ,intent(in)   ,value              :: NumVectors
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_MultiVector();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_MultiVector_Destroy ( CT_Epetra_MultiVector_ID_t * selfID );

  subroutine Epetra_MultiVector_Destroy ( selfID ) &
        bind(C,name='Epetra_MultiVector_Destroy')
    import :: FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_MultiVector_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ReplaceGlobalValue ( CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, 
  !!     int VectorIndex, double ScalarValue );

  function Epetra_MultiVector_ReplaceGlobalValue ( selfID, GlobalRow, VectorIndex, &
        ScalarValue ) result(that) bind(C,name='Epetra_MultiVector_ReplaceGlobalValue')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                  ,intent(in)   ,value              :: VectorIndex
    real(c_double)                  ,intent(in)   ,value              :: ScalarValue
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, 
  !!     double ScalarValue);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ReplaceGlobalValue_BlockPos ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue );

  function Epetra_MultiVector_ReplaceGlobalValue_BlockPos ( selfID, GlobalBlockRow, &
        BlockRowOffset, VectorIndex, ScalarValue ) result(that) &
        bind(C,name='Epetra_MultiVector_ReplaceGlobalValue_BlockPos')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: GlobalBlockRow
    integer(c_int)                  ,intent(in)   ,value              :: BlockRowOffset
    integer(c_int)                  ,intent(in)   ,value              :: VectorIndex
    real(c_double)                  ,intent(in)   ,value              :: ScalarValue
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_SumIntoGlobalValue ( CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, 
  !!     int VectorIndex, double ScalarValue );

  function Epetra_MultiVector_SumIntoGlobalValue ( selfID, GlobalRow, VectorIndex, &
        ScalarValue ) result(that) bind(C,name='Epetra_MultiVector_SumIntoGlobalValue')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                  ,intent(in)   ,value              :: VectorIndex
    real(c_double)                  ,intent(in)   ,value              :: ScalarValue
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, 
  !!     double ScalarValue);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_SumIntoGlobalValue_BlockPos ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue );

  function Epetra_MultiVector_SumIntoGlobalValue_BlockPos ( selfID, GlobalBlockRow, &
        BlockRowOffset, VectorIndex, ScalarValue ) result(that) &
        bind(C,name='Epetra_MultiVector_SumIntoGlobalValue_BlockPos')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: GlobalBlockRow
    integer(c_int)                  ,intent(in)   ,value              :: BlockRowOffset
    integer(c_int)                  ,intent(in)   ,value              :: VectorIndex
    real(c_double)                  ,intent(in)   ,value              :: ScalarValue
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceMyValue(int MyRow, int VectorIndex, double ScalarValue);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ReplaceMyValue ( CT_Epetra_MultiVector_ID_t selfID, int MyRow, 
  !!     int VectorIndex, double ScalarValue );

  function Epetra_MultiVector_ReplaceMyValue ( selfID, MyRow, VectorIndex, ScalarValue ) result(that) &
        bind(C,name='Epetra_MultiVector_ReplaceMyValue')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: MyRow
    integer(c_int)                  ,intent(in)   ,value              :: VectorIndex
    real(c_double)                  ,intent(in)   ,value              :: ScalarValue
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ReplaceMyValue_BlockPos ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue );

  function Epetra_MultiVector_ReplaceMyValue_BlockPos ( selfID, MyBlockRow, BlockRowOffset, &
        VectorIndex, ScalarValue ) result(that) &
        bind(C,name='Epetra_MultiVector_ReplaceMyValue_BlockPos')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: MyBlockRow
    integer(c_int)                  ,intent(in)   ,value              :: BlockRowOffset
    integer(c_int)                  ,intent(in)   ,value              :: VectorIndex
    real(c_double)                  ,intent(in)   ,value              :: ScalarValue
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoMyValue(int MyRow, int VectorIndex, double ScalarValue);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_SumIntoMyValue ( CT_Epetra_MultiVector_ID_t selfID, int MyRow, 
  !!     int VectorIndex, double ScalarValue );

  function Epetra_MultiVector_SumIntoMyValue ( selfID, MyRow, VectorIndex, ScalarValue ) result(that) &
        bind(C,name='Epetra_MultiVector_SumIntoMyValue')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: MyRow
    integer(c_int)                  ,intent(in)   ,value              :: VectorIndex
    real(c_double)                  ,intent(in)   ,value              :: ScalarValue
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_SumIntoMyValue_BlockPos ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue );

  function Epetra_MultiVector_SumIntoMyValue_BlockPos ( selfID, MyBlockRow, BlockRowOffset, &
        VectorIndex, ScalarValue ) result(that) &
        bind(C,name='Epetra_MultiVector_SumIntoMyValue_BlockPos')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: MyBlockRow
    integer(c_int)                  ,intent(in)   ,value              :: BlockRowOffset
    integer(c_int)                  ,intent(in)   ,value              :: VectorIndex
    real(c_double)                  ,intent(in)   ,value              :: ScalarValue
  end function


  !> <BR> Original C++ prototype:
  !! int PutScalar (double ScalarConstant);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_PutScalar ( CT_Epetra_MultiVector_ID_t selfID, double ScalarConstant );

  function Epetra_MultiVector_PutScalar ( selfID, ScalarConstant ) result(that) &
        bind(C,name='Epetra_MultiVector_PutScalar')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                  ,intent(in)   ,value              :: ScalarConstant
  end function


  !> <BR> Original C++ prototype:
  !! int Random();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Random ( CT_Epetra_MultiVector_ID_t selfID );

  function Epetra_MultiVector_Random ( selfID ) result(that) &
        bind(C,name='Epetra_MultiVector_Random')
    import :: c_int ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractCopy(double *A, int MyLDA) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ExtractCopy_Fill2DA ( CT_Epetra_MultiVector_ID_t selfID, double * A, 
  !!     int MyLDA );

  function Epetra_MultiVector_ExtractCopy_Fill2DA ( selfID, A, MyLDA ) result(that) &
        bind(C,name='Epetra_MultiVector_ExtractCopy_Fill2DA')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: A
    integer(c_int)                  ,intent(in)   ,value              :: MyLDA
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractCopy(double **ArrayOfPointers) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ExtractCopy_FillAOP ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     double ** ArrayOfPointers );

  function Epetra_MultiVector_ExtractCopy_FillAOP ( selfID, ArrayOfPointers ) result(that) &
        bind(C,name='Epetra_MultiVector_ExtractCopy_FillAOP')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: ArrayOfPointers
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractView(double **A, int *MyLDA) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ExtractView_Set2DA ( CT_Epetra_MultiVector_ID_t selfID, double ** A, 
  !!     int * MyLDA );

  function Epetra_MultiVector_ExtractView_Set2DA ( selfID, A, MyLDA ) result(that) &
        bind(C,name='Epetra_MultiVector_ExtractView_Set2DA')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: A
    integer(c_int)                                      ,dimension(*) :: MyLDA
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractView(double ***ArrayOfPointers) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ExtractView_SetAOP ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     double *** ArrayOfPointers );

  function Epetra_MultiVector_ExtractView_SetAOP ( selfID, ArrayOfPointers ) result(that) &
        bind(C,name='Epetra_MultiVector_ExtractView_SetAOP')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: ArrayOfPointers
  end function


  !> <BR> Original C++ prototype:
  !! int Dot(const Epetra_MultiVector& A, double *Result) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Dot ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     CT_Epetra_MultiVector_ID_t AID, double * Result );

  function Epetra_MultiVector_Dot ( selfID, AID, Result ) result(that) &
        bind(C,name='Epetra_MultiVector_Dot')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: AID
    real(c_double)                                      ,dimension(*) :: Result
  end function


  !> <BR> Original C++ prototype:
  !! int Abs(const Epetra_MultiVector& A);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Abs ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     CT_Epetra_MultiVector_ID_t AID );

  function Epetra_MultiVector_Abs ( selfID, AID ) result(that) &
        bind(C,name='Epetra_MultiVector_Abs')
    import :: c_int ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: AID
  end function


  !> <BR> Original C++ prototype:
  !! int Reciprocal(const Epetra_MultiVector& A);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Reciprocal ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     CT_Epetra_MultiVector_ID_t AID );

  function Epetra_MultiVector_Reciprocal ( selfID, AID ) result(that) &
        bind(C,name='Epetra_MultiVector_Reciprocal')
    import :: c_int ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: AID
  end function


  !> <BR> Original C++ prototype:
  !! int Scale(double ScalarValue);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Scale_Self ( CT_Epetra_MultiVector_ID_t selfID, double ScalarValue );

  function Epetra_MultiVector_Scale_Self ( selfID, ScalarValue ) result(that) &
        bind(C,name='Epetra_MultiVector_Scale_Self')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                  ,intent(in)   ,value              :: ScalarValue
  end function


  !> <BR> Original C++ prototype:
  !! int Scale(double ScalarA, const Epetra_MultiVector& A);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Scale ( CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  !!     CT_Epetra_MultiVector_ID_t AID );

  function Epetra_MultiVector_Scale ( selfID, ScalarA, AID ) result(that) &
        bind(C,name='Epetra_MultiVector_Scale')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                  ,intent(in)   ,value              :: ScalarA
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: AID
  end function


  !> <BR> Original C++ prototype:
  !! int Update(double ScalarA, const Epetra_MultiVector& A, double ScalarThis);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Update_WithA ( CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  !!     CT_Epetra_MultiVector_ID_t AID, double ScalarThis );

  function Epetra_MultiVector_Update_WithA ( selfID, ScalarA, AID, ScalarThis ) result(that) &
        bind(C,name='Epetra_MultiVector_Update_WithA')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                  ,intent(in)   ,value              :: ScalarA
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: AID
    real(c_double)                  ,intent(in)   ,value              :: ScalarThis
  end function


  !> <BR> Original C++ prototype:
  !! int Update(double ScalarA, const Epetra_MultiVector& A, double ScalarB, 
  !!     const Epetra_MultiVector& B, double ScalarThis);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Update_WithAB ( CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  !!     CT_Epetra_MultiVector_ID_t AID, double ScalarB, CT_Epetra_MultiVector_ID_t BID, 
  !!     double ScalarThis );

  function Epetra_MultiVector_Update_WithAB ( selfID, ScalarA, AID, ScalarB, BID, &
        ScalarThis ) result(that) bind(C,name='Epetra_MultiVector_Update_WithAB')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                  ,intent(in)   ,value              :: ScalarA
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: AID
    real(c_double)                  ,intent(in)   ,value              :: ScalarB
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
    real(c_double)                  ,intent(in)   ,value              :: ScalarThis
  end function


  !> <BR> Original C++ prototype:
  !! int Norm1 (double * Result) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Norm1 ( CT_Epetra_MultiVector_ID_t selfID, double * Result );

  function Epetra_MultiVector_Norm1 ( selfID, Result ) result(that) &
        bind(C,name='Epetra_MultiVector_Norm1')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: Result
  end function


  !> <BR> Original C++ prototype:
  !! int Norm2 (double * Result) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Norm2 ( CT_Epetra_MultiVector_ID_t selfID, double * Result );

  function Epetra_MultiVector_Norm2 ( selfID, Result ) result(that) &
        bind(C,name='Epetra_MultiVector_Norm2')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: Result
  end function


  !> <BR> Original C++ prototype:
  !! int NormInf (double * Result) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_NormInf ( CT_Epetra_MultiVector_ID_t selfID, double * Result );

  function Epetra_MultiVector_NormInf ( selfID, Result ) result(that) &
        bind(C,name='Epetra_MultiVector_NormInf')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: Result
  end function


  !> <BR> Original C++ prototype:
  !! int NormWeighted (const Epetra_MultiVector& Weights, double * Result) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_NormWeighted ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     CT_Epetra_MultiVector_ID_t WeightsID, double * Result );

  function Epetra_MultiVector_NormWeighted ( selfID, WeightsID, Result ) result(that) &
        bind(C,name='Epetra_MultiVector_NormWeighted')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: WeightsID
    real(c_double)                                      ,dimension(*) :: Result
  end function


  !> <BR> Original C++ prototype:
  !! int MinValue (double * Result) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_MinValue ( CT_Epetra_MultiVector_ID_t selfID, double * Result );

  function Epetra_MultiVector_MinValue ( selfID, Result ) result(that) &
        bind(C,name='Epetra_MultiVector_MinValue')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: Result
  end function


  !> <BR> Original C++ prototype:
  !! int MaxValue (double * Result) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_MaxValue ( CT_Epetra_MultiVector_ID_t selfID, double * Result );

  function Epetra_MultiVector_MaxValue ( selfID, Result ) result(that) &
        bind(C,name='Epetra_MultiVector_MaxValue')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: Result
  end function


  !> <BR> Original C++ prototype:
  !! int MeanValue (double * Result) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_MeanValue ( CT_Epetra_MultiVector_ID_t selfID, double * Result );

  function Epetra_MultiVector_MeanValue ( selfID, Result ) result(that) &
        bind(C,name='Epetra_MultiVector_MeanValue')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                      ,dimension(*) :: Result
  end function


  !> <BR> Original C++ prototype:
  !! int Multiply(char TransA, char TransB, double ScalarAB, const Epetra_MultiVector& A, 
  !!     const Epetra_MultiVector& B, double ScalarThis );
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Multiply_Matrix ( CT_Epetra_MultiVector_ID_t selfID, char TransA, 
  !!     char TransB, double ScalarAB, CT_Epetra_MultiVector_ID_t AID, 
  !!     CT_Epetra_MultiVector_ID_t BID, double ScalarThis );

  function Epetra_MultiVector_Multiply_Matrix ( selfID, TransA, TransB, ScalarAB, AID, BID, &
        ScalarThis ) result(that) bind(C,name='Epetra_MultiVector_Multiply_Matrix')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_char ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)          ,intent(in)   ,value              :: TransA
    character(kind=c_char)          ,intent(in)   ,value              :: TransB
    real(c_double)                  ,intent(in)   ,value              :: ScalarAB
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
    real(c_double)                  ,intent(in)   ,value              :: ScalarThis
  end function


  !> <BR> Original C++ prototype:
  !! int Multiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B, 
  !!     double ScalarThis );
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Multiply_ByEl ( CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  !!     CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, double ScalarThis );

  function Epetra_MultiVector_Multiply_ByEl ( selfID, ScalarAB, AID, BID, ScalarThis ) result(that) &
        bind(C,name='Epetra_MultiVector_Multiply_ByEl')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                  ,intent(in)   ,value              :: ScalarAB
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
    real(c_double)                  ,intent(in)   ,value              :: ScalarThis
  end function


  !> <BR> Original C++ prototype:
  !! int ReciprocalMultiply(double ScalarAB, const Epetra_MultiVector& A, 
  !!     const Epetra_MultiVector& B, double ScalarThis );
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ReciprocalMultiply ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     double ScalarAB, CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  !!     double ScalarThis );

  function Epetra_MultiVector_ReciprocalMultiply ( selfID, ScalarAB, AID, BID, ScalarThis ) result(that) &
        bind(C,name='Epetra_MultiVector_ReciprocalMultiply')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                  ,intent(in)   ,value              :: ScalarAB
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: BID
    real(c_double)                  ,intent(in)   ,value              :: ScalarThis
  end function


  !> <BR> Original C++ prototype:
  !! int SetSeed(unsigned int Seed_in);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_SetSeed ( CT_Epetra_MultiVector_ID_t selfID, unsigned int Seed_in );

  function Epetra_MultiVector_SetSeed ( selfID, Seed_in ) result(that) &
        bind(C,name='Epetra_MultiVector_SetSeed')
    import :: c_int ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: Seed_in
  end function


  !> <BR> Original C++ prototype:
  !! unsigned int Seed();
  !> <BR> <BR> CTrilinos prototype:
  !! unsigned int Epetra_MultiVector_Seed ( CT_Epetra_MultiVector_ID_t selfID );

  function Epetra_MultiVector_Seed ( selfID ) result(that) &
        bind(C,name='Epetra_MultiVector_Seed')
    import :: c_int ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector& operator = (const Epetra_MultiVector& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_MultiVector_Assign ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     CT_Epetra_MultiVector_ID_t SourceID );

  subroutine Epetra_MultiVector_Assign ( selfID, SourceID ) &
        bind(C,name='Epetra_MultiVector_Assign')
    import :: FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: SourceID
  end subroutine


  !> <BR> Original C++ prototype:
  !! double * const & operator [] (int i) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double * Epetra_MultiVector_getArray ( CT_Epetra_MultiVector_ID_t selfID, int i );

  function Epetra_MultiVector_getArray ( selfID, i ) result(that) &
        bind(C,name='Epetra_MultiVector_getArray')
    import :: c_ptr ,FT_Epetra_MultiVector_ID_t ,c_int
    
    type(c_ptr)                                                       :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: i
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Vector * & operator () (int i) const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Vector_ID_t Epetra_MultiVector_getVector ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     int i );

  function Epetra_MultiVector_getVector ( selfID, i ) result(that) &
        bind(C,name='Epetra_MultiVector_getVector')
    import :: FT_Epetra_Vector_ID_t ,FT_Epetra_MultiVector_ID_t ,c_int
    
    type(FT_Epetra_Vector_ID_t)                                       :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: i
  end function


  !> <BR> Original C++ prototype:
  !! int NumVectors() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_NumVectors ( CT_Epetra_MultiVector_ID_t selfID );

  function Epetra_MultiVector_NumVectors ( selfID ) result(that) &
        bind(C,name='Epetra_MultiVector_NumVectors')
    import :: c_int ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MyLength() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_MyLength ( CT_Epetra_MultiVector_ID_t selfID );

  function Epetra_MultiVector_MyLength ( selfID ) result(that) &
        bind(C,name='Epetra_MultiVector_MyLength')
    import :: c_int ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int GlobalLength() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_GlobalLength ( CT_Epetra_MultiVector_ID_t selfID );

  function Epetra_MultiVector_GlobalLength ( selfID ) result(that) &
        bind(C,name='Epetra_MultiVector_GlobalLength')
    import :: c_int ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int Stride() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_Stride ( CT_Epetra_MultiVector_ID_t selfID );

  function Epetra_MultiVector_Stride ( selfID ) result(that) &
        bind(C,name='Epetra_MultiVector_Stride')
    import :: c_int ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool ConstantStride() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_MultiVector_ConstantStride ( CT_Epetra_MultiVector_ID_t selfID );

  function Epetra_MultiVector_ConstantStride ( selfID ) result(that) &
        bind(C,name='Epetra_MultiVector_ConstantStride')
    import :: FT_boolean_t ,FT_Epetra_MultiVector_ID_t
    
    integer(FT_boolean_t)                                             :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceMap(const Epetra_BlockMap& map);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MultiVector_ReplaceMap ( CT_Epetra_MultiVector_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t mapID );

  function Epetra_MultiVector_ReplaceMap ( selfID, mapID ) result(that) &
        bind(C,name='Epetra_MultiVector_ReplaceMap')
    import :: c_int ,FT_Epetra_MultiVector_ID_t ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t)   ,intent(in)   ,value              :: mapID
  end function


!> @}


!> @name Epetra_OffsetIndex interface
!! @{

  ! _________________ Epetra_OffsetIndex interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_OffsetIndex_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_OffsetIndex_Degeneralize')
    import :: FT_Epetra_OffsetIndex_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_OffsetIndex_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_OffsetIndex_Generalize ( CT_Epetra_OffsetIndex_ID_t id );

  function Epetra_OffsetIndex_Generalize ( id ) result(that) &
        bind(C,name='Epetra_OffsetIndex_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_OffsetIndex_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_OffsetIndex_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph, const Epetra_CrsGraph & TargetGraph, 
  !!     Epetra_Import & Importer );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromImporter ( CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  !!     CT_Epetra_CrsGraph_ID_t TargetGraphID, CT_Epetra_Import_ID_t ImporterID );

  function Epetra_OffsetIndex_Create_FromImporter ( SourceGraphID, TargetGraphID, &
        ImporterID ) result(that) bind(C,name='Epetra_OffsetIndex_Create_FromImporter')
    import :: FT_Epetra_OffsetIndex_ID_t ,FT_Epetra_CrsGraph_ID_t ,FT_Epetra_Import_ID_t
    
    type(FT_Epetra_OffsetIndex_ID_t)                                  :: that
    type(FT_Epetra_CrsGraph_ID_t)   ,intent(in)   ,value              :: SourceGraphID
    type(FT_Epetra_CrsGraph_ID_t)   ,intent(in)   ,value              :: TargetGraphID
    type(FT_Epetra_Import_ID_t)     ,intent(in)   ,value              :: ImporterID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph, const Epetra_CrsGraph & TargetGraph, 
  !!     Epetra_Export & Exporter );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromExporter ( CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  !!     CT_Epetra_CrsGraph_ID_t TargetGraphID, CT_Epetra_Export_ID_t ExporterID );

  function Epetra_OffsetIndex_Create_FromExporter ( SourceGraphID, TargetGraphID, &
        ExporterID ) result(that) bind(C,name='Epetra_OffsetIndex_Create_FromExporter')
    import :: FT_Epetra_OffsetIndex_ID_t ,FT_Epetra_CrsGraph_ID_t ,FT_Epetra_Export_ID_t
    
    type(FT_Epetra_OffsetIndex_ID_t)                                  :: that
    type(FT_Epetra_CrsGraph_ID_t)   ,intent(in)   ,value              :: SourceGraphID
    type(FT_Epetra_CrsGraph_ID_t)   ,intent(in)   ,value              :: TargetGraphID
    type(FT_Epetra_Export_ID_t)     ,intent(in)   ,value              :: ExporterID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_OffsetIndex(const Epetra_OffsetIndex & Indexor);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Duplicate ( CT_Epetra_OffsetIndex_ID_t IndexorID );

  function Epetra_OffsetIndex_Duplicate ( IndexorID ) result(that) &
        bind(C,name='Epetra_OffsetIndex_Duplicate')
    import :: FT_Epetra_OffsetIndex_ID_t
    
    type(FT_Epetra_OffsetIndex_ID_t)                                  :: that
    type(FT_Epetra_OffsetIndex_ID_t),intent(in)   ,value              :: IndexorID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_OffsetIndex(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_OffsetIndex_Destroy ( CT_Epetra_OffsetIndex_ID_t * selfID );

  subroutine Epetra_OffsetIndex_Destroy ( selfID ) &
        bind(C,name='Epetra_OffsetIndex_Destroy')
    import :: FT_Epetra_OffsetIndex_ID_t
    
    type(FT_Epetra_OffsetIndex_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int ** SameOffsets() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int ** Epetra_OffsetIndex_SameOffsets ( CT_Epetra_OffsetIndex_ID_t selfID );

  function Epetra_OffsetIndex_SameOffsets ( selfID ) result(that) &
        bind(C,name='Epetra_OffsetIndex_SameOffsets')
    import :: c_ptr ,FT_Epetra_OffsetIndex_ID_t
    
    type(c_ptr)                                                       :: that
    type(FT_Epetra_OffsetIndex_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ** PermuteOffsets() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int ** Epetra_OffsetIndex_PermuteOffsets ( CT_Epetra_OffsetIndex_ID_t selfID );

  function Epetra_OffsetIndex_PermuteOffsets ( selfID ) result(that) &
        bind(C,name='Epetra_OffsetIndex_PermuteOffsets')
    import :: c_ptr ,FT_Epetra_OffsetIndex_ID_t
    
    type(c_ptr)                                                       :: that
    type(FT_Epetra_OffsetIndex_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ** RemoteOffsets() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int ** Epetra_OffsetIndex_RemoteOffsets ( CT_Epetra_OffsetIndex_ID_t selfID );

  function Epetra_OffsetIndex_RemoteOffsets ( selfID ) result(that) &
        bind(C,name='Epetra_OffsetIndex_RemoteOffsets')
    import :: c_ptr ,FT_Epetra_OffsetIndex_ID_t
    
    type(c_ptr)                                                       :: that
    type(FT_Epetra_OffsetIndex_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_Object interface
!! @{

  ! _________________ Epetra_Object interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Object_ID_t Epetra_Object_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Object_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Object_Degeneralize')
    import :: FT_Epetra_Object_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Object_ID_t)                                   :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Object_Generalize ( CT_Epetra_Object_ID_t id );

  function Epetra_Object_Generalize ( id ) result(that) &
        bind(C,name='Epetra_Object_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Object_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Object_ID_t) ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Object(int TracebackModeIn = -1, bool set_label = true);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Object_ID_t Epetra_Object_Create ( int TracebackModeIn, boolean set_label );

  function Epetra_Object_Create ( TracebackModeIn, set_label ) result(that) &
        bind(C,name='Epetra_Object_Create')
    import :: FT_Epetra_Object_ID_t ,c_int ,FT_boolean_t
    
    type(FT_Epetra_Object_ID_t)                                   :: that
    integer(c_int)              ,intent(in)   ,value              :: TracebackModeIn
    integer(FT_boolean_t)       ,intent(in)   ,value              :: set_label
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Object(const char * const Label, int TracebackModeIn = -1);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Object_ID_t Epetra_Object_Create_WithLabel ( const char * const Label, 
  !!     int TracebackModeIn );

  function Epetra_Object_Create_WithLabel ( Label, TracebackModeIn ) result(that) &
        bind(C,name='Epetra_Object_Create_WithLabel')
    import :: FT_Epetra_Object_ID_t ,c_char ,c_int
    
    type(FT_Epetra_Object_ID_t)                                   :: that
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: Label
    integer(c_int)              ,intent(in)   ,value              :: TracebackModeIn
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Object(const Epetra_Object& Object);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Object_ID_t Epetra_Object_Duplicate ( CT_Epetra_Object_ID_t ObjectID );

  function Epetra_Object_Duplicate ( ObjectID ) result(that) &
        bind(C,name='Epetra_Object_Duplicate')
    import :: FT_Epetra_Object_ID_t
    
    type(FT_Epetra_Object_ID_t)                                   :: that
    type(FT_Epetra_Object_ID_t) ,intent(in)   ,value              :: ObjectID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Object();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Object_Destroy ( CT_Epetra_Object_ID_t * selfID );

  subroutine Epetra_Object_Destroy ( selfID ) bind(C,name='Epetra_Object_Destroy')
    import :: FT_Epetra_Object_ID_t
    
    type(FT_Epetra_Object_ID_t)                                   :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual void SetLabel(const char * const Label);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Object_SetLabel ( CT_Epetra_Object_ID_t selfID, const char * const Label );

  subroutine Epetra_Object_SetLabel ( selfID, Label ) bind(C,name='Epetra_Object_SetLabel')
    import :: FT_Epetra_Object_ID_t ,c_char
    
    type(FT_Epetra_Object_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: Label
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual const char * Label() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Epetra_Object_Label ( CT_Epetra_Object_ID_t selfID );

  function Epetra_Object_Label ( selfID ) result(that) bind(C,name='Epetra_Object_Label')
    import :: c_ptr ,FT_Epetra_Object_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Object_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! static void SetTracebackMode(int TracebackModeValue);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Object_SetTracebackMode ( int TracebackModeValue );

  subroutine Epetra_Object_SetTracebackMode ( TracebackModeValue ) &
        bind(C,name='Epetra_Object_SetTracebackMode')
    import :: c_int
    
    integer(c_int)              ,intent(in)   ,value              :: TracebackModeValue
  end subroutine


  !> <BR> Original C++ prototype:
  !! static int GetTracebackMode();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Object_GetTracebackMode (  );

  function Epetra_Object_GetTracebackMode (  ) result(that) &
        bind(C,name='Epetra_Object_GetTracebackMode')
    import :: c_int
    
    integer(c_int)                                                :: that
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ReportError(const string Message, int ErrorCode) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Object_ReportError ( CT_Epetra_Object_ID_t selfID, const char Message[], 
  !!     int ErrorCode );

  function Epetra_Object_ReportError ( selfID, Message, ErrorCode ) result(that) &
        bind(C,name='Epetra_Object_ReportError')
    import :: c_int ,FT_Epetra_Object_ID_t ,c_char
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Object_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: Message
    integer(c_int)              ,intent(in)   ,value              :: ErrorCode
  end function


!> @}


!> @name Epetra_RowMatrix interface
!! @{

  ! _________________ Epetra_RowMatrix interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_RowMatrix_ID_t Epetra_RowMatrix_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_RowMatrix_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_RowMatrix_Degeneralize')
    import :: FT_Epetra_RowMatrix_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_RowMatrix_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_RowMatrix_Generalize ( CT_Epetra_RowMatrix_ID_t id );

  function Epetra_RowMatrix_Generalize ( id ) result(that) &
        bind(C,name='Epetra_RowMatrix_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_RowMatrix_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_RowMatrix();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_RowMatrix_Destroy ( CT_Epetra_RowMatrix_ID_t * selfID );

  subroutine Epetra_RowMatrix_Destroy ( selfID ) bind(C,name='Epetra_RowMatrix_Destroy')
    import :: FT_Epetra_RowMatrix_ID_t
    
    type(FT_Epetra_RowMatrix_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual int NumMyRowEntries(int MyRow, int & NumEntries) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_NumMyRowEntries ( CT_Epetra_RowMatrix_ID_t selfID, int MyRow, 
  !!     int * NumEntries );

  function Epetra_RowMatrix_NumMyRowEntries ( selfID, MyRow, NumEntries ) result(that) &
        bind(C,name='Epetra_RowMatrix_NumMyRowEntries')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(inout)                    :: NumEntries
  end function


  !> <BR> Original C++ prototype:
  !! virtual int MaxNumEntries() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_MaxNumEntries ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_MaxNumEntries ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_MaxNumEntries')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, 
  !!     int * Indices) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_ExtractMyRowCopy ( CT_Epetra_RowMatrix_ID_t selfID, int MyRow, 
  !!     int Length, int * NumEntries, double * Values, int * Indices );

  function Epetra_RowMatrix_ExtractMyRowCopy ( selfID, MyRow, Length, NumEntries, Values, &
        Indices ) result(that) bind(C,name='Epetra_RowMatrix_ExtractMyRowCopy')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(in)   ,value              :: Length
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_ExtractDiagonalCopy ( CT_Epetra_RowMatrix_ID_t selfID, 
  !!     CT_Epetra_Vector_ID_t DiagonalID );

  function Epetra_RowMatrix_ExtractDiagonalCopy ( selfID, DiagonalID ) result(that) &
        bind(C,name='Epetra_RowMatrix_ExtractDiagonalCopy')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: DiagonalID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 
  !!     0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_Multiply ( CT_Epetra_RowMatrix_ID_t selfID, boolean TransA, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

  function Epetra_RowMatrix_Multiply ( selfID, TransA, XID, YID ) result(that) &
        bind(C,name='Epetra_RowMatrix_Multiply')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t ,FT_boolean_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: TransA
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
  !!     Epetra_MultiVector& Y) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_Solve ( CT_Epetra_RowMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  !!     boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

  function Epetra_RowMatrix_Solve ( selfID, Upper, Trans, UnitDiagonal, XID, YID ) result(that) &
        bind(C,name='Epetra_RowMatrix_Solve')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t ,FT_boolean_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: Upper
    integer(FT_boolean_t)         ,intent(in)   ,value              :: Trans
    integer(FT_boolean_t)         ,intent(in)   ,value              :: UnitDiagonal
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int InvRowSums(Epetra_Vector& x) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_InvRowSums ( CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_RowMatrix_InvRowSums ( selfID, xID ) result(that) &
        bind(C,name='Epetra_RowMatrix_InvRowSums')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int LeftScale(const Epetra_Vector& x) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_LeftScale ( CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_RowMatrix_LeftScale ( selfID, xID ) result(that) &
        bind(C,name='Epetra_RowMatrix_LeftScale')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int InvColSums(Epetra_Vector& x) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_InvColSums ( CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_RowMatrix_InvColSums ( selfID, xID ) result(that) &
        bind(C,name='Epetra_RowMatrix_InvColSums')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int RightScale(const Epetra_Vector& x) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_RightScale ( CT_Epetra_RowMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_RowMatrix_RightScale ( selfID, xID ) result(that) &
        bind(C,name='Epetra_RowMatrix_RightScale')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool Filled() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_RowMatrix_Filled ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_Filled ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_Filled')
    import :: FT_boolean_t ,FT_Epetra_RowMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double NormInf() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_RowMatrix_NormInf ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NormInf ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NormInf')
    import :: c_double ,FT_Epetra_RowMatrix_ID_t
    
    real(c_double)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double NormOne() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_RowMatrix_NormOne ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NormOne ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NormOne')
    import :: c_double ,FT_Epetra_RowMatrix_ID_t
    
    real(c_double)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumGlobalNonzeros() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_NumGlobalNonzeros ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NumGlobalNonzeros ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NumGlobalNonzeros')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumGlobalRows() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_NumGlobalRows ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NumGlobalRows ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NumGlobalRows')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumGlobalCols() const= 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_NumGlobalCols ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NumGlobalCols ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NumGlobalCols')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumGlobalDiagonals() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_NumGlobalDiagonals ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NumGlobalDiagonals ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NumGlobalDiagonals')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumMyNonzeros() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_NumMyNonzeros ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NumMyNonzeros ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NumMyNonzeros')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumMyRows() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_NumMyRows ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NumMyRows ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NumMyRows')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumMyCols() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_NumMyCols ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NumMyCols ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NumMyCols')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumMyDiagonals() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_RowMatrix_NumMyDiagonals ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_NumMyDiagonals ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_NumMyDiagonals')
    import :: c_int ,FT_Epetra_RowMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool LowerTriangular() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_RowMatrix_LowerTriangular ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_LowerTriangular ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_LowerTriangular')
    import :: FT_boolean_t ,FT_Epetra_RowMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool UpperTriangular() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_RowMatrix_UpperTriangular ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_UpperTriangular ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_UpperTriangular')
    import :: FT_boolean_t ,FT_Epetra_RowMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_Map & RowMatrixRowMap() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixRowMap ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_RowMatrixRowMap ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_RowMatrixRowMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_RowMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_Map & RowMatrixColMap() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_RowMatrix_RowMatrixColMap ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_RowMatrixColMap ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_RowMatrixColMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_RowMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_Import * RowMatrixImporter() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Import_ID_t Epetra_RowMatrix_RowMatrixImporter ( CT_Epetra_RowMatrix_ID_t selfID );

  function Epetra_RowMatrix_RowMatrixImporter ( selfID ) result(that) &
        bind(C,name='Epetra_RowMatrix_RowMatrixImporter')
    import :: FT_Epetra_Import_ID_t ,FT_Epetra_RowMatrix_ID_t
    
    type(FT_Epetra_Import_ID_t)                                     :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_CompObject interface
!! @{

  ! _________________ Epetra_CompObject interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_CompObject_ID_t Epetra_CompObject_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_CompObject_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_CompObject_Degeneralize')
    import :: FT_Epetra_CompObject_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_CompObject_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_CompObject_Generalize ( CT_Epetra_CompObject_ID_t id );

  function Epetra_CompObject_Generalize ( id ) result(that) &
        bind(C,name='Epetra_CompObject_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_CompObject_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CompObject();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CompObject_ID_t Epetra_CompObject_Create (  );

  function Epetra_CompObject_Create (  ) result(that) &
        bind(C,name='Epetra_CompObject_Create')
    import :: FT_Epetra_CompObject_ID_t
    
    type(FT_Epetra_CompObject_ID_t)                                  :: that
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CompObject(const Epetra_CompObject& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CompObject_ID_t Epetra_CompObject_Duplicate ( CT_Epetra_CompObject_ID_t SourceID );

  function Epetra_CompObject_Duplicate ( SourceID ) result(that) &
        bind(C,name='Epetra_CompObject_Duplicate')
    import :: FT_Epetra_CompObject_ID_t
    
    type(FT_Epetra_CompObject_ID_t)                                  :: that
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: SourceID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_CompObject();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_Destroy ( CT_Epetra_CompObject_ID_t * selfID );

  subroutine Epetra_CompObject_Destroy ( selfID ) bind(C,name='Epetra_CompObject_Destroy')
    import :: FT_Epetra_CompObject_ID_t
    
    type(FT_Epetra_CompObject_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SetFlopCounter(const Epetra_Flops & FlopCounter_in);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_SetFlopCounter ( CT_Epetra_CompObject_ID_t selfID, 
  !!     CT_Epetra_Flops_ID_t FlopCounter_inID );

  subroutine Epetra_CompObject_SetFlopCounter ( selfID, FlopCounter_inID ) &
        bind(C,name='Epetra_CompObject_SetFlopCounter')
    import :: FT_Epetra_CompObject_ID_t ,FT_Epetra_Flops_ID_t
    
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Flops_ID_t)     ,intent(in)   ,value              :: FlopCounter_inID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SetFlopCounter(const Epetra_CompObject & CompObject);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_SetFlopCounter_Matching ( CT_Epetra_CompObject_ID_t selfID, 
  !!     CT_Epetra_CompObject_ID_t CompObjectID );

  subroutine Epetra_CompObject_SetFlopCounter_Matching ( selfID, CompObjectID ) &
        bind(C,name='Epetra_CompObject_SetFlopCounter_Matching')
    import :: FT_Epetra_CompObject_ID_t
    
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: CompObjectID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void UnsetFlopCounter();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_UnsetFlopCounter ( CT_Epetra_CompObject_ID_t selfID );

  subroutine Epetra_CompObject_UnsetFlopCounter ( selfID ) &
        bind(C,name='Epetra_CompObject_UnsetFlopCounter')
    import :: FT_Epetra_CompObject_ID_t
    
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! Epetra_Flops * GetFlopCounter() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Flops_ID_t Epetra_CompObject_GetFlopCounter ( CT_Epetra_CompObject_ID_t selfID );

  function Epetra_CompObject_GetFlopCounter ( selfID ) result(that) &
        bind(C,name='Epetra_CompObject_GetFlopCounter')
    import :: FT_Epetra_Flops_ID_t ,FT_Epetra_CompObject_ID_t
    
    type(FT_Epetra_Flops_ID_t)                                       :: that
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void ResetFlops() const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_ResetFlops ( CT_Epetra_CompObject_ID_t selfID );

  subroutine Epetra_CompObject_ResetFlops ( selfID ) &
        bind(C,name='Epetra_CompObject_ResetFlops')
    import :: FT_Epetra_CompObject_ID_t
    
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! double Flops() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_CompObject_Flops ( CT_Epetra_CompObject_ID_t selfID );

  function Epetra_CompObject_Flops ( selfID ) result(that) &
        bind(C,name='Epetra_CompObject_Flops')
    import :: c_double ,FT_Epetra_CompObject_ID_t
    
    real(c_double)                                                   :: that
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void UpdateFlops(int Flops_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_UpdateFlops_Int ( CT_Epetra_CompObject_ID_t selfID, int Flops_in );

  subroutine Epetra_CompObject_UpdateFlops_Int ( selfID, Flops_in ) &
        bind(C,name='Epetra_CompObject_UpdateFlops_Int')
    import :: FT_Epetra_CompObject_ID_t ,c_int
    
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                 ,intent(in)   ,value              :: Flops_in
  end subroutine


  !> <BR> Original C++ prototype:
  !! void UpdateFlops(long int Flops_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_UpdateFlops_Long ( CT_Epetra_CompObject_ID_t selfID, 
  !!     long int Flops_in );

  subroutine Epetra_CompObject_UpdateFlops_Long ( selfID, Flops_in ) &
        bind(C,name='Epetra_CompObject_UpdateFlops_Long')
    import :: FT_Epetra_CompObject_ID_t ,c_int
    
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                 ,intent(in)   ,value              :: Flops_in
  end subroutine


  !> <BR> Original C++ prototype:
  !! void UpdateFlops(double Flops_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_UpdateFlops_Double ( CT_Epetra_CompObject_ID_t selfID, 
  !!     double Flops_in );

  subroutine Epetra_CompObject_UpdateFlops_Double ( selfID, Flops_in ) &
        bind(C,name='Epetra_CompObject_UpdateFlops_Double')
    import :: FT_Epetra_CompObject_ID_t ,c_double
    
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                 ,intent(in)   ,value              :: Flops_in
  end subroutine


  !> <BR> Original C++ prototype:
  !! void UpdateFlops(float Flops_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_UpdateFlops_Float ( CT_Epetra_CompObject_ID_t selfID, float Flops_in );

  subroutine Epetra_CompObject_UpdateFlops_Float ( selfID, Flops_in ) &
        bind(C,name='Epetra_CompObject_UpdateFlops_Float')
    import :: FT_Epetra_CompObject_ID_t ,c_float
    
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
    real(c_float)                  ,intent(in)   ,value              :: Flops_in
  end subroutine


  !> <BR> Original C++ prototype:
  !! Epetra_CompObject& operator=(const Epetra_CompObject& src);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CompObject_Assign ( CT_Epetra_CompObject_ID_t selfID, 
  !!     CT_Epetra_CompObject_ID_t srcID );

  subroutine Epetra_CompObject_Assign ( selfID, srcID ) &
        bind(C,name='Epetra_CompObject_Assign')
    import :: FT_Epetra_CompObject_ID_t
    
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_CompObject_ID_t),intent(in)   ,value              :: srcID
  end subroutine


!> @}


!> @name Epetra_Directory interface
!! @{

  ! _________________ Epetra_Directory interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Directory_ID_t Epetra_Directory_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Directory_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Directory_Degeneralize')
    import :: FT_Epetra_Directory_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Directory_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Directory_Generalize ( CT_Epetra_Directory_ID_t id );

  function Epetra_Directory_Generalize ( id ) result(that) &
        bind(C,name='Epetra_Directory_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Directory_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Directory_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Directory();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Directory_Destroy ( CT_Epetra_Directory_ID_t * selfID );

  subroutine Epetra_Directory_Destroy ( selfID ) bind(C,name='Epetra_Directory_Destroy')
    import :: FT_Epetra_Directory_ID_t
    
    type(FT_Epetra_Directory_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual int GetDirectoryEntries( const Epetra_BlockMap& Map, const int NumEntries, 
  !!     const int * GlobalEntries, int * Procs, int * LocalEntries, int * EntrySizes, 
  !!     bool high_rank_sharing_procs=false) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Directory_GetDirectoryEntries ( CT_Epetra_Directory_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t MapID, const int NumEntries, const int * GlobalEntries, 
  !!     int * Procs, int * LocalEntries, int * EntrySizes, boolean high_rank_sharing_procs );

  function Epetra_Directory_GetDirectoryEntries ( selfID, MapID, NumEntries, GlobalEntries, &
        Procs, LocalEntries, EntrySizes, high_rank_sharing_procs ) result(that) &
        bind(C,name='Epetra_Directory_GetDirectoryEntries')
    import :: c_int ,FT_Epetra_Directory_ID_t ,FT_Epetra_BlockMap_ID_t ,FT_boolean_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_Directory_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t) ,intent(in)   ,value              :: MapID
    integer(c_int)                ,intent(in)   ,value              :: NumEntries
    integer(c_int)                ,intent(in)         ,dimension(*) :: GlobalEntries
    integer(c_int)                                    ,dimension(*) :: Procs
    integer(c_int)                                    ,dimension(*) :: LocalEntries
    integer(c_int)                                    ,dimension(*) :: EntrySizes
    integer(FT_boolean_t)         ,intent(in)   ,value              :: high_rank_sharing_procs
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool GIDsAllUniquelyOwned() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_Directory_GIDsAllUniquelyOwned ( CT_Epetra_Directory_ID_t selfID );

  function Epetra_Directory_GIDsAllUniquelyOwned ( selfID ) result(that) &
        bind(C,name='Epetra_Directory_GIDsAllUniquelyOwned')
    import :: FT_boolean_t ,FT_Epetra_Directory_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_Directory_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_Flops interface
!! @{

  ! _________________ Epetra_Flops interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Flops_ID_t Epetra_Flops_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Flops_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Flops_Degeneralize')
    import :: FT_Epetra_Flops_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Flops_ID_t)                                    :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Flops_Generalize ( CT_Epetra_Flops_ID_t id );

  function Epetra_Flops_Generalize ( id ) result(that) &
        bind(C,name='Epetra_Flops_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Flops_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Flops_ID_t)  ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Flops(void);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Flops_ID_t Epetra_Flops_Create (  );

  function Epetra_Flops_Create (  ) result(that) bind(C,name='Epetra_Flops_Create')
    import :: FT_Epetra_Flops_ID_t
    
    type(FT_Epetra_Flops_ID_t)                                    :: that
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Flops(const Epetra_Flops& Flops_in);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( CT_Epetra_Flops_ID_t Flops_inID );

  function Epetra_Flops_Duplicate ( Flops_inID ) result(that) &
        bind(C,name='Epetra_Flops_Duplicate')
    import :: FT_Epetra_Flops_ID_t
    
    type(FT_Epetra_Flops_ID_t)                                    :: that
    type(FT_Epetra_Flops_ID_t)  ,intent(in)   ,value              :: Flops_inID
  end function


  !> <BR> Original C++ prototype:
  !! double Flops() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID );

  function Epetra_Flops_Flops ( selfID ) result(that) bind(C,name='Epetra_Flops_Flops')
    import :: c_double ,FT_Epetra_Flops_ID_t
    
    real(c_double)                                                :: that
    type(FT_Epetra_Flops_ID_t)  ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void ResetFlops();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID );

  subroutine Epetra_Flops_ResetFlops ( selfID ) bind(C,name='Epetra_Flops_ResetFlops')
    import :: FT_Epetra_Flops_ID_t
    
    type(FT_Epetra_Flops_ID_t)  ,intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Flops(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID );

  subroutine Epetra_Flops_Destroy ( selfID ) bind(C,name='Epetra_Flops_Destroy')
    import :: FT_Epetra_Flops_ID_t
    
    type(FT_Epetra_Flops_ID_t)                                    :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! Epetra_Flops& operator=(const Epetra_Flops& src);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Flops_Assign ( CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID );

  subroutine Epetra_Flops_Assign ( selfID, srcID ) bind(C,name='Epetra_Flops_Assign')
    import :: FT_Epetra_Flops_ID_t
    
    type(FT_Epetra_Flops_ID_t)  ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Flops_ID_t)  ,intent(in)   ,value              :: srcID
  end subroutine


!> @}


!> @name Epetra_SrcDistObject interface
!! @{

  ! _________________ Epetra_SrcDistObject interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_SrcDistObject_ID_t Epetra_SrcDistObject_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_SrcDistObject_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_SrcDistObject_Degeneralize')
    import :: FT_Epetra_SrcDistObject_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_SrcDistObject_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t)  ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_SrcDistObject_Generalize ( CT_Epetra_SrcDistObject_ID_t id );

  function Epetra_SrcDistObject_Generalize ( id ) result(that) &
        bind(C,name='Epetra_SrcDistObject_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_SrcDistObject_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                    :: that
    type(FT_Epetra_SrcDistObject_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_SrcDistObject();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SrcDistObject_Destroy ( CT_Epetra_SrcDistObject_ID_t * selfID );

  subroutine Epetra_SrcDistObject_Destroy ( selfID ) &
        bind(C,name='Epetra_SrcDistObject_Destroy')
    import :: FT_Epetra_SrcDistObject_ID_t
    
    type(FT_Epetra_SrcDistObject_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_BlockMap & Map() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_SrcDistObject_Map ( CT_Epetra_SrcDistObject_ID_t selfID );

  function Epetra_SrcDistObject_Map ( selfID ) result(that) &
        bind(C,name='Epetra_SrcDistObject_Map')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_SrcDistObject_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                       :: that
    type(FT_Epetra_SrcDistObject_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_MpiComm interface
!! @{

  ! _________________ Epetra_MpiComm interface bodies _________________





#ifdef HAVE_MPI


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_MpiComm_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_MpiComm_Degeneralize')
    import :: FT_Epetra_MpiComm_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_MpiComm_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_MpiComm_Generalize ( CT_Epetra_MpiComm_ID_t id );

  function Epetra_MpiComm_Generalize ( id ) result(that) &
        bind(C,name='Epetra_MpiComm_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_MpiComm_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MpiComm(const Epetra_MpiComm & Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( CT_Epetra_MpiComm_ID_t CommID );

  function Epetra_MpiComm_Duplicate ( CommID ) result(that) &
        bind(C,name='Epetra_MpiComm_Duplicate')
    import :: FT_Epetra_MpiComm_ID_t
    
    type(FT_Epetra_MpiComm_ID_t)                                  :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Comm * Clone() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Epetra_MpiComm_Clone ( CT_Epetra_MpiComm_ID_t selfID );

  function Epetra_MpiComm_Clone ( selfID ) result(that) bind(C,name='Epetra_MpiComm_Clone')
    import :: FT_Epetra_Comm_ID_t ,FT_Epetra_MpiComm_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                     :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_MpiComm();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID );

  subroutine Epetra_MpiComm_Destroy ( selfID ) bind(C,name='Epetra_MpiComm_Destroy')
    import :: FT_Epetra_MpiComm_ID_t
    
    type(FT_Epetra_MpiComm_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void Barrier() const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_MpiComm_Barrier ( CT_Epetra_MpiComm_ID_t selfID );

  subroutine Epetra_MpiComm_Barrier ( selfID ) bind(C,name='Epetra_MpiComm_Barrier')
    import :: FT_Epetra_MpiComm_ID_t
    
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int Broadcast(double * MyVals, int Count, int Root) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_Broadcast_Double ( CT_Epetra_MpiComm_ID_t selfID, double * MyVals, 
  !!     int Count, int Root );

  function Epetra_MpiComm_Broadcast_Double ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_MpiComm_Broadcast_Double')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: MyVals
    integer(c_int)              ,intent(in)   ,value              :: Count
    integer(c_int)              ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! int Broadcast(int * MyVals, int Count, int Root) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_Broadcast_Int ( CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int Count, 
  !!     int Root );

  function Epetra_MpiComm_Broadcast_Int ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_MpiComm_Broadcast_Int')
    import :: c_int ,FT_Epetra_MpiComm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: MyVals
    integer(c_int)              ,intent(in)   ,value              :: Count
    integer(c_int)              ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! int Broadcast(long * MyVals, int Count, int Root) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_Broadcast_Long ( CT_Epetra_MpiComm_ID_t selfID, long * MyVals, int Count, 
  !!     int Root );

  function Epetra_MpiComm_Broadcast_Long ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_MpiComm_Broadcast_Long')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: MyVals
    integer(c_int)              ,intent(in)   ,value              :: Count
    integer(c_int)              ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! int Broadcast(char * MyVals, int Count, int Root) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_Broadcast_Char ( CT_Epetra_MpiComm_ID_t selfID, char * MyVals, int Count, 
  !!     int Root );

  function Epetra_MpiComm_Broadcast_Char ( selfID, MyVals, Count, Root ) result(that) &
        bind(C,name='Epetra_MpiComm_Broadcast_Char')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_char
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                          ,dimension(*) :: MyVals
    integer(c_int)              ,intent(in)   ,value              :: Count
    integer(c_int)              ,intent(in)   ,value              :: Root
  end function


  !> <BR> Original C++ prototype:
  !! int GatherAll(double * MyVals, double * AllVals, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_GatherAll_Double ( CT_Epetra_MpiComm_ID_t selfID, double * MyVals, 
  !!     double * AllVals, int Count );

  function Epetra_MpiComm_GatherAll_Double ( selfID, MyVals, AllVals, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_GatherAll_Double')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: MyVals
    real(c_double)                                  ,dimension(*) :: AllVals
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int GatherAll(int * MyVals, int * AllVals, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_GatherAll_Int ( CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * AllVals, 
  !!     int Count );

  function Epetra_MpiComm_GatherAll_Int ( selfID, MyVals, AllVals, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_GatherAll_Int')
    import :: c_int ,FT_Epetra_MpiComm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: MyVals
    integer(c_int)                                  ,dimension(*) :: AllVals
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int GatherAll(long * MyVals, long * AllVals, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_GatherAll_Long ( CT_Epetra_MpiComm_ID_t selfID, long * MyVals, 
  !!     long * AllVals, int Count );

  function Epetra_MpiComm_GatherAll_Long ( selfID, MyVals, AllVals, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_GatherAll_Long')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: MyVals
    integer(c_long)                                 ,dimension(*) :: AllVals
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int SumAll(double * PartialSums, double * GlobalSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_SumAll_Double ( CT_Epetra_MpiComm_ID_t selfID, double * PartialSums, 
  !!     double * GlobalSums, int Count );

  function Epetra_MpiComm_SumAll_Double ( selfID, PartialSums, GlobalSums, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_SumAll_Double')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: PartialSums
    real(c_double)                                  ,dimension(*) :: GlobalSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int SumAll(int * PartialSums, int * GlobalSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_SumAll_Int ( CT_Epetra_MpiComm_ID_t selfID, int * PartialSums, 
  !!     int * GlobalSums, int Count );

  function Epetra_MpiComm_SumAll_Int ( selfID, PartialSums, GlobalSums, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_SumAll_Int')
    import :: c_int ,FT_Epetra_MpiComm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: PartialSums
    integer(c_int)                                  ,dimension(*) :: GlobalSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int SumAll(long * PartialSums, long * GlobalSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_SumAll_Long ( CT_Epetra_MpiComm_ID_t selfID, long * PartialSums, 
  !!     long * GlobalSums, int Count );

  function Epetra_MpiComm_SumAll_Long ( selfID, PartialSums, GlobalSums, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_SumAll_Long')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: PartialSums
    integer(c_long)                                 ,dimension(*) :: GlobalSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_MaxAll_Double ( CT_Epetra_MpiComm_ID_t selfID, double * PartialMaxs, 
  !!     double * GlobalMaxs, int Count );

  function Epetra_MpiComm_MaxAll_Double ( selfID, PartialMaxs, GlobalMaxs, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_MaxAll_Double')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: PartialMaxs
    real(c_double)                                  ,dimension(*) :: GlobalMaxs
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_MaxAll_Int ( CT_Epetra_MpiComm_ID_t selfID, int * PartialMaxs, 
  !!     int * GlobalMaxs, int Count );

  function Epetra_MpiComm_MaxAll_Int ( selfID, PartialMaxs, GlobalMaxs, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_MaxAll_Int')
    import :: c_int ,FT_Epetra_MpiComm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: PartialMaxs
    integer(c_int)                                  ,dimension(*) :: GlobalMaxs
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_MaxAll_Long ( CT_Epetra_MpiComm_ID_t selfID, long * PartialMaxs, 
  !!     long * GlobalMaxs, int Count );

  function Epetra_MpiComm_MaxAll_Long ( selfID, PartialMaxs, GlobalMaxs, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_MaxAll_Long')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: PartialMaxs
    integer(c_long)                                 ,dimension(*) :: GlobalMaxs
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MinAll(double * PartialMins, double * GlobalMins, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_MinAll_Double ( CT_Epetra_MpiComm_ID_t selfID, double * PartialMins, 
  !!     double * GlobalMins, int Count );

  function Epetra_MpiComm_MinAll_Double ( selfID, PartialMins, GlobalMins, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_MinAll_Double')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: PartialMins
    real(c_double)                                  ,dimension(*) :: GlobalMins
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MinAll(int * PartialMins, int * GlobalMins, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_MinAll_Int ( CT_Epetra_MpiComm_ID_t selfID, int * PartialMins, 
  !!     int * GlobalMins, int Count );

  function Epetra_MpiComm_MinAll_Int ( selfID, PartialMins, GlobalMins, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_MinAll_Int')
    import :: c_int ,FT_Epetra_MpiComm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: PartialMins
    integer(c_int)                                  ,dimension(*) :: GlobalMins
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MinAll(long * PartialMins, long * GlobalMins, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_MinAll_Long ( CT_Epetra_MpiComm_ID_t selfID, long * PartialMins, 
  !!     long * GlobalMins, int Count );

  function Epetra_MpiComm_MinAll_Long ( selfID, PartialMins, GlobalMins, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_MinAll_Long')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: PartialMins
    integer(c_long)                                 ,dimension(*) :: GlobalMins
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int ScanSum(double * MyVals, double * ScanSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_ScanSum_Double ( CT_Epetra_MpiComm_ID_t selfID, double * MyVals, 
  !!     double * ScanSums, int Count );

  function Epetra_MpiComm_ScanSum_Double ( selfID, MyVals, ScanSums, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_ScanSum_Double')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: MyVals
    real(c_double)                                  ,dimension(*) :: ScanSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int ScanSum(int * MyVals, int * ScanSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_ScanSum_Int ( CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * ScanSums, 
  !!     int Count );

  function Epetra_MpiComm_ScanSum_Int ( selfID, MyVals, ScanSums, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_ScanSum_Int')
    import :: c_int ,FT_Epetra_MpiComm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                  ,dimension(*) :: MyVals
    integer(c_int)                                  ,dimension(*) :: ScanSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int ScanSum(long * MyVals, long * ScanSums, int Count) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_ScanSum_Long ( CT_Epetra_MpiComm_ID_t selfID, long * MyVals, 
  !!     long * ScanSums, int Count );

  function Epetra_MpiComm_ScanSum_Long ( selfID, MyVals, ScanSums, Count ) result(that) &
        bind(C,name='Epetra_MpiComm_ScanSum_Long')
    import :: c_int ,FT_Epetra_MpiComm_ID_t ,c_long
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    integer(c_long)                                 ,dimension(*) :: MyVals
    integer(c_long)                                 ,dimension(*) :: ScanSums
    integer(c_int)              ,intent(in)   ,value              :: Count
  end function


  !> <BR> Original C++ prototype:
  !! int MyPID() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_MyPID ( CT_Epetra_MpiComm_ID_t selfID );

  function Epetra_MpiComm_MyPID ( selfID ) result(that) bind(C,name='Epetra_MpiComm_MyPID')
    import :: c_int ,FT_Epetra_MpiComm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumProc() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_NumProc ( CT_Epetra_MpiComm_ID_t selfID );

  function Epetra_MpiComm_NumProc ( selfID ) result(that) &
        bind(C,name='Epetra_MpiComm_NumProc')
    import :: c_int ,FT_Epetra_MpiComm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Distributor * CreateDistributor() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Distributor_ID_t Epetra_MpiComm_CreateDistributor ( CT_Epetra_MpiComm_ID_t selfID );

  function Epetra_MpiComm_CreateDistributor ( selfID ) result(that) &
        bind(C,name='Epetra_MpiComm_CreateDistributor')
    import :: FT_Epetra_Distributor_ID_t ,FT_Epetra_MpiComm_ID_t
    
    type(FT_Epetra_Distributor_ID_t)                                  :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Directory * CreateDirectory(const Epetra_BlockMap & Map) const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Directory_ID_t Epetra_MpiComm_CreateDirectory ( CT_Epetra_MpiComm_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t MapID );

  function Epetra_MpiComm_CreateDirectory ( selfID, MapID ) result(that) &
        bind(C,name='Epetra_MpiComm_CreateDirectory')
    import :: FT_Epetra_Directory_ID_t ,FT_Epetra_MpiComm_ID_t ,FT_Epetra_BlockMap_ID_t
    
    type(FT_Epetra_Directory_ID_t)                                  :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: MapID
  end function


  !> <BR> Original C++ prototype:
  !! int GetMpiTag() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_MpiComm_GetMpiTag ( CT_Epetra_MpiComm_ID_t selfID );

  function Epetra_MpiComm_GetMpiTag ( selfID ) result(that) &
        bind(C,name='Epetra_MpiComm_GetMpiTag')
    import :: c_int ,FT_Epetra_MpiComm_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MpiComm & operator=(const Epetra_MpiComm & Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_MpiComm_Assign ( CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_MpiComm_ID_t CommID );

  subroutine Epetra_MpiComm_Assign ( selfID, CommID ) bind(C,name='Epetra_MpiComm_Assign')
    import :: FT_Epetra_MpiComm_ID_t
    
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MpiComm_ID_t),intent(in)   ,value              :: CommID
  end subroutine


#endif /* HAVE_MPI */


!> @}


!> @name Epetra_CrsMatrix interface
!! @{

  ! _________________ Epetra_CrsMatrix interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_CrsMatrix_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Degeneralize')
    import :: FT_Epetra_CrsMatrix_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_CrsMatrix_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_CrsMatrix_Generalize ( CT_Epetra_CrsMatrix_ID_t id );

  function Epetra_CrsMatrix_Generalize ( id ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const int* NumEntriesPerRow, 
  !!     bool StaticProfile = false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_Map_ID_t RowMapID, const int * NumEntriesPerRow, boolean StaticProfile );

  function Epetra_CrsMatrix_Create_VarPerRow ( CV, RowMapID, NumEntriesPerRow, &
        StaticProfile ) result(that) bind(C,name='Epetra_CrsMatrix_Create_VarPerRow')
    import :: FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_Map_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_CrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: RowMapID
    integer(c_int)                ,intent(in)         ,dimension(*) :: NumEntriesPerRow
    integer(FT_boolean_t)         ,intent(in)   ,value              :: StaticProfile
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int NumEntriesPerRow, 
  !!     bool StaticProfile = false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_Map_ID_t RowMapID, int NumEntriesPerRow, boolean StaticProfile );

  function Epetra_CrsMatrix_Create ( CV, RowMapID, NumEntriesPerRow, StaticProfile ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Create')
    import :: FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_Map_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_CrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: RowMapID
    integer(c_int)                ,intent(in)   ,value              :: NumEntriesPerRow
    integer(FT_boolean_t)         ,intent(in)   ,value              :: StaticProfile
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, 
  !!     const int* NumEntriesPerRow, bool StaticProfile = false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_Map_ID_t RowMapID, CT_Epetra_Map_ID_t ColMapID, const int * NumEntriesPerRow, 
  !!     boolean StaticProfile );

  function Epetra_CrsMatrix_Create_VarPerRow_WithColMap ( CV, RowMapID, ColMapID, &
        NumEntriesPerRow, StaticProfile ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Create_VarPerRow_WithColMap')
    import :: FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_Map_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_CrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: RowMapID
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: ColMapID
    integer(c_int)                ,intent(in)         ,dimension(*) :: NumEntriesPerRow
    integer(FT_boolean_t)         ,intent(in)   ,value              :: StaticProfile
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, 
  !!     int NumEntriesPerRow, bool StaticProfile = false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_WithColMap ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_Map_ID_t RowMapID, CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  !!     boolean StaticProfile );

  function Epetra_CrsMatrix_Create_WithColMap ( CV, RowMapID, ColMapID, NumEntriesPerRow, &
        StaticProfile ) result(that) bind(C,name='Epetra_CrsMatrix_Create_WithColMap')
    import :: FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_Map_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_CrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: RowMapID
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: ColMapID
    integer(c_int)                ,intent(in)   ,value              :: NumEntriesPerRow
    integer(FT_boolean_t)         ,intent(in)   ,value              :: StaticProfile
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& Graph);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_FromGraph ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_CrsGraph_ID_t GraphID );

  function Epetra_CrsMatrix_Create_FromGraph ( CV, GraphID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Create_FromGraph')
    import :: FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_CrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_CrsGraph_ID_t) ,intent(in)   ,value              :: GraphID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsMatrix(const Epetra_CrsMatrix& Matrix);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( CT_Epetra_CrsMatrix_ID_t MatrixID );

  function Epetra_CrsMatrix_Duplicate ( MatrixID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Duplicate')
    import :: FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_CrsMatrix_ID_t)                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: MatrixID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_CrsMatrix();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CrsMatrix_Destroy ( CT_Epetra_CrsMatrix_ID_t * selfID );

  subroutine Epetra_CrsMatrix_Destroy ( selfID ) bind(C,name='Epetra_CrsMatrix_Destroy')
    import :: FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_CrsMatrix_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! Epetra_CrsMatrix& operator=(const Epetra_CrsMatrix& src);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CrsMatrix_Assign ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     CT_Epetra_CrsMatrix_ID_t srcID );

  subroutine Epetra_CrsMatrix_Assign ( selfID, srcID ) &
        bind(C,name='Epetra_CrsMatrix_Assign')
    import :: FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: srcID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int PutScalar(double ScalarConstant);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_PutScalar ( CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant );

  function Epetra_CrsMatrix_PutScalar ( selfID, ScalarConstant ) result(that) &
        bind(C,name='Epetra_CrsMatrix_PutScalar')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                ,intent(in)   ,value              :: ScalarConstant
  end function


  !> <BR> Original C++ prototype:
  !! int Scale(double ScalarConstant);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_Scale ( CT_Epetra_CrsMatrix_ID_t selfID, double ScalarConstant );

  function Epetra_CrsMatrix_Scale ( selfID, ScalarConstant ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Scale')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                ,intent(in)   ,value              :: ScalarConstant
  end function


  !> <BR> Original C++ prototype:
  !! virtual int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_InsertGlobalValues ( CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, 
  !!     int NumEntries, double * Values, int * Indices );

  function Epetra_CrsMatrix_InsertGlobalValues ( selfID, GlobalRow, NumEntries, Values, &
        Indices ) result(that) bind(C,name='Epetra_CrsMatrix_InsertGlobalValues')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                ,intent(in)   ,value              :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ReplaceGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ReplaceGlobalValues ( CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, 
  !!     int NumEntries, double * Values, int * Indices );

  function Epetra_CrsMatrix_ReplaceGlobalValues ( selfID, GlobalRow, NumEntries, Values, &
        Indices ) result(that) bind(C,name='Epetra_CrsMatrix_ReplaceGlobalValues')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                ,intent(in)   ,value              :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! virtual int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_SumIntoGlobalValues ( CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, 
  !!     int NumEntries, double * Values, int * Indices );

  function Epetra_CrsMatrix_SumIntoGlobalValues ( selfID, GlobalRow, NumEntries, Values, &
        Indices ) result(that) bind(C,name='Epetra_CrsMatrix_SumIntoGlobalValues')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                ,intent(in)   ,value              :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int InsertMyValues(int MyRow, int NumEntries, double* Values, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_InsertMyValues ( CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, 
  !!     int NumEntries, double * Values, int * Indices );

  function Epetra_CrsMatrix_InsertMyValues ( selfID, MyRow, NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_CrsMatrix_InsertMyValues')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(in)   ,value              :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceMyValues(int MyRow, int NumEntries, double* Values, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ReplaceMyValues ( CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, 
  !!     int NumEntries, double * Values, int * Indices );

  function Epetra_CrsMatrix_ReplaceMyValues ( selfID, MyRow, NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ReplaceMyValues')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(in)   ,value              :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoMyValues(int MyRow, int NumEntries, double* Values, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_SumIntoMyValues ( CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, 
  !!     int NumEntries, double * Values, int * Indices );

  function Epetra_CrsMatrix_SumIntoMyValues ( selfID, MyRow, NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_CrsMatrix_SumIntoMyValues')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(in)   ,value              :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceDiagonalValues(const Epetra_Vector& Diagonal);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ReplaceDiagonalValues ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     CT_Epetra_Vector_ID_t DiagonalID );

  function Epetra_CrsMatrix_ReplaceDiagonalValues ( selfID, DiagonalID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ReplaceDiagonalValues')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: DiagonalID
  end function


  !> <BR> Original C++ prototype:
  !! int FillComplete(bool OptimizeDataStorage = true);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_FillComplete ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     boolean OptimizeDataStorage );

  function Epetra_CrsMatrix_FillComplete ( selfID, OptimizeDataStorage ) result(that) &
        bind(C,name='Epetra_CrsMatrix_FillComplete')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_boolean_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: OptimizeDataStorage
  end function


  !> <BR> Original C++ prototype:
  !! int FillComplete(const Epetra_Map& DomainMap, const Epetra_Map& RangeMap, 
  !!     bool OptimizeDataStorage = true);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_FillComplete_UsingMaps ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     CT_Epetra_Map_ID_t DomainMapID, CT_Epetra_Map_ID_t RangeMapID, 
  !!     boolean OptimizeDataStorage );

  function Epetra_CrsMatrix_FillComplete_UsingMaps ( selfID, DomainMapID, RangeMapID, &
        OptimizeDataStorage ) result(that) &
        bind(C,name='Epetra_CrsMatrix_FillComplete_UsingMaps')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Map_ID_t ,FT_boolean_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: DomainMapID
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: RangeMapID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: OptimizeDataStorage
  end function


  !> <BR> Original C++ prototype:
  !! int OptimizeStorage();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_OptimizeStorage ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_OptimizeStorage ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_OptimizeStorage')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MakeDataContiguous();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_MakeDataContiguous ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_MakeDataContiguous ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_MakeDataContiguous')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values, 
  !!     int* Indices) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     int GlobalRow, int Length, int * NumEntries, double * Values, int * Indices );

  function Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices ( selfID, GlobalRow, Length, &
        NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                ,intent(in)   ,value              :: Length
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values, 
  !!     int* Indices) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     int MyRow, int Length, int * NumEntries, double * Values, int * Indices );

  function Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices ( selfID, MyRow, Length, &
        NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ExtractMyRowCopy_WithIndices')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(in)   ,value              :: Length
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractGlobalRowCopy(int GlobalRow, int Length, int& NumEntries, double* Values) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ExtractGlobalRowCopy ( CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, 
  !!     int Length, int * NumEntries, double * Values );

  function Epetra_CrsMatrix_ExtractGlobalRowCopy ( selfID, GlobalRow, Length, NumEntries, &
        Values ) result(that) bind(C,name='Epetra_CrsMatrix_ExtractGlobalRowCopy')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                ,intent(in)   ,value              :: Length
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, double* Values) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ExtractMyRowCopy ( CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, 
  !!     int Length, int * NumEntries, double * Values );

  function Epetra_CrsMatrix_ExtractMyRowCopy ( selfID, MyRow, Length, NumEntries, Values ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ExtractMyRowCopy')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(in)   ,value              :: Length
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractDiagonalCopy(Epetra_Vector& Diagonal) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ExtractDiagonalCopy ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     CT_Epetra_Vector_ID_t DiagonalID );

  function Epetra_CrsMatrix_ExtractDiagonalCopy ( selfID, DiagonalID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ExtractDiagonalCopy')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: DiagonalID
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values, int*& Indices) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     int GlobalRow, int * NumEntries, double ** Values, int ** Indices );

  function Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices ( selfID, GlobalRow, &
        NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ExtractGlobalRowView_WithIndices')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                ,intent(inout)      ,dimension(*) :: Values
    integer(c_int)                ,intent(inout)      ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractMyRowView(int MyRow, int& NumEntries, double*& Values, int*& Indices) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ExtractMyRowView_WithIndices ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     int MyRow, int * NumEntries, double ** Values, int ** Indices );

  function Epetra_CrsMatrix_ExtractMyRowView_WithIndices ( selfID, MyRow, NumEntries, &
        Values, Indices ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ExtractMyRowView_WithIndices')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                ,intent(inout)      ,dimension(*) :: Values
    integer(c_int)                ,intent(inout)      ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractGlobalRowView(int GlobalRow, int& NumEntries, double*& Values) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ExtractGlobalRowView ( CT_Epetra_CrsMatrix_ID_t selfID, int GlobalRow, 
  !!     int * NumEntries, double ** Values );

  function Epetra_CrsMatrix_ExtractGlobalRowView ( selfID, GlobalRow, NumEntries, Values ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ExtractGlobalRowView')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                ,intent(inout)      ,dimension(*) :: Values
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractMyRowView(int MyRow, int& NumEntries, double*& Values) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ExtractMyRowView ( CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, 
  !!     int * NumEntries, double ** Values );

  function Epetra_CrsMatrix_ExtractMyRowView ( selfID, MyRow, NumEntries, Values ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ExtractMyRowView')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                ,intent(inout)      ,dimension(*) :: Values
  end function


  !> <BR> Original C++ prototype:
  !! int Multiply(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_Multiply_Vector ( CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  !!     CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID );

  function Epetra_CrsMatrix_Multiply_Vector ( selfID, TransA, xID, yID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Multiply_Vector')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_boolean_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: TransA
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: yID
  end function


  !> <BR> Original C++ prototype:
  !! int Multiply1(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_Multiply1_Vector ( CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  !!     CT_Epetra_Vector_ID_t xID, CT_Epetra_Vector_ID_t yID );

  function Epetra_CrsMatrix_Multiply1_Vector ( selfID, TransA, xID, yID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Multiply1_Vector')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_boolean_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: TransA
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: yID
  end function


  !> <BR> Original C++ prototype:
  !! int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_Multiply_MultiVector ( CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

  function Epetra_CrsMatrix_Multiply_MultiVector ( selfID, TransA, XID, YID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Multiply_MultiVector')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_boolean_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: TransA
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! int Multiply1(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_Multiply1_MultiVector ( CT_Epetra_CrsMatrix_ID_t selfID, boolean TransA, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

  function Epetra_CrsMatrix_Multiply1_MultiVector ( selfID, TransA, XID, YID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Multiply1_MultiVector')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_boolean_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: TransA
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, 
  !!     Epetra_Vector& y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_Solve_Vector ( CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, 
  !!     boolean Trans, boolean UnitDiagonal, CT_Epetra_Vector_ID_t xID, 
  !!     CT_Epetra_Vector_ID_t yID );

  function Epetra_CrsMatrix_Solve_Vector ( selfID, Upper, Trans, UnitDiagonal, xID, yID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Solve_Vector')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_boolean_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: Upper
    integer(FT_boolean_t)         ,intent(in)   ,value              :: Trans
    integer(FT_boolean_t)         ,intent(in)   ,value              :: UnitDiagonal
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: yID
  end function


  !> <BR> Original C++ prototype:
  !! int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
  !!     Epetra_MultiVector& Y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_Solve_MultiVector ( CT_Epetra_CrsMatrix_ID_t selfID, boolean Upper, 
  !!     boolean Trans, boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, 
  !!     CT_Epetra_MultiVector_ID_t YID );

  function Epetra_CrsMatrix_Solve_MultiVector ( selfID, Upper, Trans, UnitDiagonal, XID, &
        YID ) result(that) bind(C,name='Epetra_CrsMatrix_Solve_MultiVector')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_boolean_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: Upper
    integer(FT_boolean_t)         ,intent(in)   ,value              :: Trans
    integer(FT_boolean_t)         ,intent(in)   ,value              :: UnitDiagonal
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! int InvRowSums(Epetra_Vector& x) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_InvRowSums ( CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_CrsMatrix_InvRowSums ( selfID, xID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_InvRowSums')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! int InvRowMaxs(Epetra_Vector& x) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_InvRowMaxs ( CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_CrsMatrix_InvRowMaxs ( selfID, xID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_InvRowMaxs')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! int LeftScale(const Epetra_Vector& x);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_LeftScale ( CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_CrsMatrix_LeftScale ( selfID, xID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_LeftScale')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! int InvColSums(Epetra_Vector& x) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_InvColSums ( CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_CrsMatrix_InvColSums ( selfID, xID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_InvColSums')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! int InvColMaxs(Epetra_Vector& x) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_InvColMaxs ( CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_CrsMatrix_InvColMaxs ( selfID, xID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_InvColMaxs')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! int RightScale(const Epetra_Vector& x);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_RightScale ( CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_Vector_ID_t xID );

  function Epetra_CrsMatrix_RightScale ( selfID, xID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_RightScale')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)   ,intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! bool Filled() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_Filled ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_Filled ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Filled')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool StorageOptimized() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_StorageOptimized ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_StorageOptimized ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_StorageOptimized')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool IndicesAreGlobal() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_IndicesAreGlobal ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_IndicesAreGlobal ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_IndicesAreGlobal')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool IndicesAreLocal() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_IndicesAreLocal ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_IndicesAreLocal ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_IndicesAreLocal')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool IndicesAreContiguous() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_IndicesAreContiguous ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_IndicesAreContiguous ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_IndicesAreContiguous')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool LowerTriangular() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_LowerTriangular ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_LowerTriangular ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_LowerTriangular')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool UpperTriangular() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_UpperTriangular ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_UpperTriangular ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_UpperTriangular')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool NoDiagonal() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_NoDiagonal ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NoDiagonal ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NoDiagonal')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double NormInf() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_CrsMatrix_NormInf ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NormInf ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NormInf')
    import :: c_double ,FT_Epetra_CrsMatrix_ID_t
    
    real(c_double)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double NormOne() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_CrsMatrix_NormOne ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NormOne ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NormOne')
    import :: c_double ,FT_Epetra_CrsMatrix_ID_t
    
    real(c_double)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double NormFrobenius() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_CrsMatrix_NormFrobenius ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NormFrobenius ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NormFrobenius')
    import :: c_double ,FT_Epetra_CrsMatrix_ID_t
    
    real(c_double)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalNonzeros() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumGlobalNonzeros ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NumGlobalNonzeros ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumGlobalNonzeros')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalRows() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumGlobalRows ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NumGlobalRows ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumGlobalRows')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalCols() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumGlobalCols ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NumGlobalCols ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumGlobalCols')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalDiagonals() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumGlobalDiagonals ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NumGlobalDiagonals ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumGlobalDiagonals')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyNonzeros() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumMyNonzeros ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NumMyNonzeros ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumMyNonzeros')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyRows() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumMyRows ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NumMyRows ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumMyRows')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyCols() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumMyCols ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NumMyCols ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumMyCols')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyDiagonals() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumMyDiagonals ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_NumMyDiagonals ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumMyDiagonals')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalEntries(int Row) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumGlobalEntries ( CT_Epetra_CrsMatrix_ID_t selfID, int Row );

  function Epetra_CrsMatrix_NumGlobalEntries ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumGlobalEntries')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int NumAllocatedGlobalEntries(int Row) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumAllocatedGlobalEntries ( CT_Epetra_CrsMatrix_ID_t selfID, int Row );

  function Epetra_CrsMatrix_NumAllocatedGlobalEntries ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumAllocatedGlobalEntries')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int MaxNumEntries() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_MaxNumEntries ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_MaxNumEntries ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_MaxNumEntries')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int GlobalMaxNumEntries() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_GlobalMaxNumEntries ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_GlobalMaxNumEntries ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_GlobalMaxNumEntries')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyEntries(int Row) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumMyEntries ( CT_Epetra_CrsMatrix_ID_t selfID, int Row );

  function Epetra_CrsMatrix_NumMyEntries ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumMyEntries')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int NumAllocatedMyEntries(int Row) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumAllocatedMyEntries ( CT_Epetra_CrsMatrix_ID_t selfID, int Row );

  function Epetra_CrsMatrix_NumAllocatedMyEntries ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumAllocatedMyEntries')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int IndexBase() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_IndexBase ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_IndexBase ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_IndexBase')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool StaticGraph();
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_StaticGraph ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_StaticGraph ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_StaticGraph')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_CrsGraph& Graph() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsGraph_ID_t Epetra_CrsMatrix_Graph ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_Graph ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Graph')
    import :: FT_Epetra_CrsGraph_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_CrsGraph_ID_t)                                   :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Map& RowMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_RowMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_RowMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceRowMap(const Epetra_BlockMap& newmap);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ReplaceRowMap ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t newmapID );

  function Epetra_CrsMatrix_ReplaceRowMap ( selfID, newmapID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ReplaceRowMap')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t) ,intent(in)   ,value              :: newmapID
  end function


  !> <BR> Original C++ prototype:
  !! bool HaveColMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_HaveColMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_HaveColMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_HaveColMap')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceColMap(const Epetra_BlockMap& newmap);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ReplaceColMap ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t newmapID );

  function Epetra_CrsMatrix_ReplaceColMap ( selfID, newmapID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ReplaceColMap')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t) ,intent(in)   ,value              :: newmapID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Map& ColMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_CrsMatrix_ColMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_ColMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ColMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Map& DomainMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_CrsMatrix_DomainMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_DomainMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_DomainMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Map& RangeMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_CrsMatrix_RangeMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_RangeMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_RangeMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Import* Importer() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Import_ID_t Epetra_CrsMatrix_Importer ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_Importer ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Importer')
    import :: FT_Epetra_Import_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Import_ID_t)                                     :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Export* Exporter() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Export_ID_t Epetra_CrsMatrix_Exporter ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_Exporter ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Exporter')
    import :: FT_Epetra_Export_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Export_ID_t)                                     :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Comm& Comm() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Epetra_CrsMatrix_Comm ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_Comm ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Comm')
    import :: FT_Epetra_Comm_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                       :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int LRID( int GRID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_LRID ( CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in );

  function Epetra_CrsMatrix_LRID ( selfID, GRID_in ) result(that) &
        bind(C,name='Epetra_CrsMatrix_LRID')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GRID_in
  end function


  !> <BR> Original C++ prototype:
  !! int GRID( int LRID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_GRID ( CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in );

  function Epetra_CrsMatrix_GRID ( selfID, LRID_in ) result(that) &
        bind(C,name='Epetra_CrsMatrix_GRID')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: LRID_in
  end function


  !> <BR> Original C++ prototype:
  !! int LCID( int GCID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_LCID ( CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in );

  function Epetra_CrsMatrix_LCID ( selfID, GCID_in ) result(that) &
        bind(C,name='Epetra_CrsMatrix_LCID')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GCID_in
  end function


  !> <BR> Original C++ prototype:
  !! int GCID( int LCID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_GCID ( CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in );

  function Epetra_CrsMatrix_GCID ( selfID, LCID_in ) result(that) &
        bind(C,name='Epetra_CrsMatrix_GCID')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: LCID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyGRID(int GRID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_MyGRID ( CT_Epetra_CrsMatrix_ID_t selfID, int GRID_in );

  function Epetra_CrsMatrix_MyGRID ( selfID, GRID_in ) result(that) &
        bind(C,name='Epetra_CrsMatrix_MyGRID')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t ,c_int
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GRID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyLRID(int LRID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_MyLRID ( CT_Epetra_CrsMatrix_ID_t selfID, int LRID_in );

  function Epetra_CrsMatrix_MyLRID ( selfID, LRID_in ) result(that) &
        bind(C,name='Epetra_CrsMatrix_MyLRID')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t ,c_int
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: LRID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyGCID(int GCID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_MyGCID ( CT_Epetra_CrsMatrix_ID_t selfID, int GCID_in );

  function Epetra_CrsMatrix_MyGCID ( selfID, GCID_in ) result(that) &
        bind(C,name='Epetra_CrsMatrix_MyGCID')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t ,c_int
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GCID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyLCID(int LCID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_MyLCID ( CT_Epetra_CrsMatrix_ID_t selfID, int LCID_in );

  function Epetra_CrsMatrix_MyLCID ( selfID, LCID_in ) result(that) &
        bind(C,name='Epetra_CrsMatrix_MyLCID')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t ,c_int
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: LCID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyGlobalRow(int GID) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_MyGlobalRow ( CT_Epetra_CrsMatrix_ID_t selfID, int GID );

  function Epetra_CrsMatrix_MyGlobalRow ( selfID, GID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_MyGlobalRow')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t ,c_int
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: GID
  end function


  !> <BR> Original C++ prototype:
  !! const char* Label() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Epetra_CrsMatrix_Label ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_Label ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Label')
    import :: c_ptr ,FT_Epetra_CrsMatrix_ID_t
    
    type(c_ptr)                                                     :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int SetUseTranspose(bool UseTranspose_in);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_SetUseTranspose ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     boolean UseTranspose_in );

  function Epetra_CrsMatrix_SetUseTranspose ( selfID, UseTranspose_in ) result(that) &
        bind(C,name='Epetra_CrsMatrix_SetUseTranspose')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_boolean_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: UseTranspose_in
  end function


  !> <BR> Original C++ prototype:
  !! int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_Apply ( CT_Epetra_CrsMatrix_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, 
  !!     CT_Epetra_MultiVector_ID_t YID );

  function Epetra_CrsMatrix_Apply ( selfID, XID, YID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_Apply')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_ApplyInverse ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

  function Epetra_CrsMatrix_ApplyInverse ( selfID, XID, YID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ApplyInverse')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! bool HasNormInf() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_HasNormInf ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_HasNormInf ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_HasNormInf')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool UseTranspose() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsMatrix_UseTranspose ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_UseTranspose ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_UseTranspose')
    import :: FT_boolean_t ,FT_Epetra_CrsMatrix_ID_t
    
    integer(FT_boolean_t)                                           :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Map& OperatorDomainMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorDomainMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_OperatorDomainMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_OperatorDomainMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Map& OperatorRangeMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_CrsMatrix_OperatorRangeMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_OperatorRangeMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_OperatorRangeMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyRowEntries(int MyRow, int& NumEntries) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_NumMyRowEntries ( CT_Epetra_CrsMatrix_ID_t selfID, int MyRow, 
  !!     int * NumEntries );

  function Epetra_CrsMatrix_NumMyRowEntries ( selfID, MyRow, NumEntries ) result(that) &
        bind(C,name='Epetra_CrsMatrix_NumMyRowEntries')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(inout)                    :: NumEntries
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Map& RowMatrixRowMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixRowMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_RowMatrixRowMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_RowMatrixRowMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Map& RowMatrixColMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_CrsMatrix_RowMatrixColMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_RowMatrixColMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_RowMatrixColMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Import* RowMatrixImporter() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Import_ID_t Epetra_CrsMatrix_RowMatrixImporter ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_RowMatrixImporter ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_RowMatrixImporter')
    import :: FT_Epetra_Import_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Import_ID_t)                                     :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! inline double* operator[] (int Loc) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double * Epetra_CrsMatrix_getRow ( CT_Epetra_CrsMatrix_ID_t selfID, int Loc );

  function Epetra_CrsMatrix_getRow ( selfID, Loc ) result(that) &
        bind(C,name='Epetra_CrsMatrix_getRow')
    import :: c_ptr ,FT_Epetra_CrsMatrix_ID_t ,c_int
    
    type(c_ptr)                                                     :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: Loc
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Map& ImportMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_CrsMatrix_ImportMap ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_ImportMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_ImportMap')
    import :: FT_Epetra_Map_ID_t ,FT_Epetra_CrsMatrix_ID_t
    
    type(FT_Epetra_Map_ID_t)                                        :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int TransformToLocal();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_TransformToLocal ( CT_Epetra_CrsMatrix_ID_t selfID );

  function Epetra_CrsMatrix_TransformToLocal ( selfID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_TransformToLocal')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int TransformToLocal(const Epetra_Map* DomainMap, const Epetra_Map* RangeMap);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsMatrix_TransformToLocal_UsingMaps ( CT_Epetra_CrsMatrix_ID_t selfID, 
  !!     CT_Epetra_Map_ID_t DomainMapID, CT_Epetra_Map_ID_t RangeMapID );

  function Epetra_CrsMatrix_TransformToLocal_UsingMaps ( selfID, DomainMapID, RangeMapID ) result(that) &
        bind(C,name='Epetra_CrsMatrix_TransformToLocal_UsingMaps')
    import :: c_int ,FT_Epetra_CrsMatrix_ID_t ,FT_Epetra_Map_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_CrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: DomainMapID
    type(FT_Epetra_Map_ID_t)      ,intent(in)   ,value              :: RangeMapID
  end function


!> @}


!> @name Epetra_CrsGraph interface
!! @{

  ! _________________ Epetra_CrsGraph interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_CrsGraph_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_CrsGraph_Degeneralize')
    import :: FT_Epetra_CrsGraph_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_CrsGraph_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_CrsGraph_Generalize ( CT_Epetra_CrsGraph_ID_t id );

  function Epetra_CrsGraph_Generalize ( id ) result(that) &
        bind(C,name='Epetra_CrsGraph_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_CrsGraph_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, 
  !!     const int* NumIndicesPerRow, bool StaticProfile = false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_BlockMap_ID_t RowMapID, const int * NumIndicesPerRow, boolean StaticProfile );

  function Epetra_CrsGraph_Create_VarPerRow ( CV, RowMapID, NumIndicesPerRow, StaticProfile ) result(that) &
        bind(C,name='Epetra_CrsGraph_Create_VarPerRow')
    import :: FT_Epetra_CrsGraph_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_BlockMap_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_CrsGraph_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: RowMapID
    integer(c_int)               ,intent(in)         ,dimension(*) :: NumIndicesPerRow
    integer(FT_boolean_t)        ,intent(in)   ,value              :: StaticProfile
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, int NumIndicesPerRow, 
  !!     bool StaticProfile = false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_BlockMap_ID_t RowMapID, int NumIndicesPerRow, boolean StaticProfile );

  function Epetra_CrsGraph_Create ( CV, RowMapID, NumIndicesPerRow, StaticProfile ) result(that) &
        bind(C,name='Epetra_CrsGraph_Create')
    import :: FT_Epetra_CrsGraph_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_BlockMap_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_CrsGraph_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: RowMapID
    integer(c_int)               ,intent(in)   ,value              :: NumIndicesPerRow
    integer(FT_boolean_t)        ,intent(in)   ,value              :: StaticProfile
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, 
  !!     const Epetra_BlockMap& ColMap, const int* NumIndicesPerRow, bool StaticProfile = false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_VarPerRow_WithColMap ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_BlockMap_ID_t RowMapID, CT_Epetra_BlockMap_ID_t ColMapID, 
  !!     const int * NumIndicesPerRow, boolean StaticProfile );

  function Epetra_CrsGraph_Create_VarPerRow_WithColMap ( CV, RowMapID, ColMapID, &
        NumIndicesPerRow, StaticProfile ) result(that) &
        bind(C,name='Epetra_CrsGraph_Create_VarPerRow_WithColMap')
    import :: FT_Epetra_CrsGraph_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_BlockMap_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_CrsGraph_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: RowMapID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: ColMapID
    integer(c_int)               ,intent(in)         ,dimension(*) :: NumIndicesPerRow
    integer(FT_boolean_t)        ,intent(in)   ,value              :: StaticProfile
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsGraph(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, 
  !!     const Epetra_BlockMap& ColMap, int NumIndicesPerRow, bool StaticProfile = false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Create_With_ColMap ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_BlockMap_ID_t RowMapID, CT_Epetra_BlockMap_ID_t ColMapID, int NumIndicesPerRow, 
  !!     boolean StaticProfile );

  function Epetra_CrsGraph_Create_With_ColMap ( CV, RowMapID, ColMapID, NumIndicesPerRow, &
        StaticProfile ) result(that) bind(C,name='Epetra_CrsGraph_Create_With_ColMap')
    import :: FT_Epetra_CrsGraph_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_BlockMap_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_CrsGraph_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: RowMapID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: ColMapID
    integer(c_int)               ,intent(in)   ,value              :: NumIndicesPerRow
    integer(FT_boolean_t)        ,intent(in)   ,value              :: StaticProfile
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsGraph(const Epetra_CrsGraph& Graph);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_CrsGraph_ID_t Epetra_CrsGraph_Duplicate ( CT_Epetra_CrsGraph_ID_t GraphID );

  function Epetra_CrsGraph_Duplicate ( GraphID ) result(that) &
        bind(C,name='Epetra_CrsGraph_Duplicate')
    import :: FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_CrsGraph_ID_t)                                  :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: GraphID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_CrsGraph();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CrsGraph_Destroy ( CT_Epetra_CrsGraph_ID_t * selfID );

  subroutine Epetra_CrsGraph_Destroy ( selfID ) bind(C,name='Epetra_CrsGraph_Destroy')
    import :: FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_CrsGraph_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int InsertGlobalIndices(int GlobalRow, int NumIndices, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_InsertGlobalIndices ( CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, 
  !!     int NumIndices, int * Indices );

  function Epetra_CrsGraph_InsertGlobalIndices ( selfID, GlobalRow, NumIndices, Indices ) result(that) &
        bind(C,name='Epetra_CrsGraph_InsertGlobalIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GlobalRow
    integer(c_int)               ,intent(in)   ,value              :: NumIndices
    integer(c_int)                                   ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int RemoveGlobalIndices(int GlobalRow, int NumIndices, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_RemoveGlobalIndices ( CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, 
  !!     int NumIndices, int * Indices );

  function Epetra_CrsGraph_RemoveGlobalIndices ( selfID, GlobalRow, NumIndices, Indices ) result(that) &
        bind(C,name='Epetra_CrsGraph_RemoveGlobalIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GlobalRow
    integer(c_int)               ,intent(in)   ,value              :: NumIndices
    integer(c_int)                                   ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int RemoveGlobalIndices(int Row);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_RemoveGlobalIndices_LocalRow ( CT_Epetra_CrsGraph_ID_t selfID, int Row );

  function Epetra_CrsGraph_RemoveGlobalIndices_LocalRow ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsGraph_RemoveGlobalIndices_LocalRow')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int InsertMyIndices(int LocalRow, int NumIndices, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_InsertMyIndices ( CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, 
  !!     int NumIndices, int * Indices );

  function Epetra_CrsGraph_InsertMyIndices ( selfID, LocalRow, NumIndices, Indices ) result(that) &
        bind(C,name='Epetra_CrsGraph_InsertMyIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LocalRow
    integer(c_int)               ,intent(in)   ,value              :: NumIndices
    integer(c_int)                                   ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int RemoveMyIndices(int LocalRow, int NumIndices, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_RemoveMyIndices ( CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, 
  !!     int NumIndices, int * Indices );

  function Epetra_CrsGraph_RemoveMyIndices ( selfID, LocalRow, NumIndices, Indices ) result(that) &
        bind(C,name='Epetra_CrsGraph_RemoveMyIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LocalRow
    integer(c_int)               ,intent(in)   ,value              :: NumIndices
    integer(c_int)                                   ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int RemoveMyIndices(int Row);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_RemoveMyIndices_LocalRow ( CT_Epetra_CrsGraph_ID_t selfID, int Row );

  function Epetra_CrsGraph_RemoveMyIndices_LocalRow ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsGraph_RemoveMyIndices_LocalRow')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int FillComplete();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_FillComplete ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_FillComplete ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_FillComplete')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int FillComplete(const Epetra_BlockMap& DomainMap, const Epetra_BlockMap& RangeMap);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_FillComplete_UsingMaps ( CT_Epetra_CrsGraph_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t DomainMapID, CT_Epetra_BlockMap_ID_t RangeMapID );

  function Epetra_CrsGraph_FillComplete_UsingMaps ( selfID, DomainMapID, RangeMapID ) result(that) &
        bind(C,name='Epetra_CrsGraph_FillComplete_UsingMaps')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: DomainMapID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: RangeMapID
  end function


  !> <BR> Original C++ prototype:
  !! int OptimizeStorage();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_OptimizeStorage ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_OptimizeStorage ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_OptimizeStorage')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractGlobalRowCopy(int GlobalRow, int LenOfIndices, int& NumIndices, int* Indices) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_ExtractGlobalRowCopy ( CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, 
  !!     int LenOfIndices, int * NumIndices, int * Indices );

  function Epetra_CrsGraph_ExtractGlobalRowCopy ( selfID, GlobalRow, LenOfIndices, &
        NumIndices, Indices ) result(that) &
        bind(C,name='Epetra_CrsGraph_ExtractGlobalRowCopy')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GlobalRow
    integer(c_int)               ,intent(in)   ,value              :: LenOfIndices
    integer(c_int)               ,intent(inout)                    :: NumIndices
    integer(c_int)                                   ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractMyRowCopy(int LocalRow, int LenOfIndices, int& NumIndices, int* Indices) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_ExtractMyRowCopy ( CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, 
  !!     int LenOfIndices, int * NumIndices, int * Indices );

  function Epetra_CrsGraph_ExtractMyRowCopy ( selfID, LocalRow, LenOfIndices, NumIndices, &
        Indices ) result(that) bind(C,name='Epetra_CrsGraph_ExtractMyRowCopy')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LocalRow
    integer(c_int)               ,intent(in)   ,value              :: LenOfIndices
    integer(c_int)               ,intent(inout)                    :: NumIndices
    integer(c_int)                                   ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractGlobalRowView(int GlobalRow, int& NumIndices, int*& Indices) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_ExtractGlobalRowView ( CT_Epetra_CrsGraph_ID_t selfID, int GlobalRow, 
  !!     int * NumIndices, int ** Indices );

  function Epetra_CrsGraph_ExtractGlobalRowView ( selfID, GlobalRow, NumIndices, Indices ) result(that) &
        bind(C,name='Epetra_CrsGraph_ExtractGlobalRowView')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GlobalRow
    integer(c_int)               ,intent(inout)                    :: NumIndices
    integer(c_int)               ,intent(inout)      ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractMyRowView(int LocalRow, int& NumIndices, int*& Indices) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_ExtractMyRowView ( CT_Epetra_CrsGraph_ID_t selfID, int LocalRow, 
  !!     int * NumIndices, int ** Indices );

  function Epetra_CrsGraph_ExtractMyRowView ( selfID, LocalRow, NumIndices, Indices ) result(that) &
        bind(C,name='Epetra_CrsGraph_ExtractMyRowView')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LocalRow
    integer(c_int)               ,intent(inout)                    :: NumIndices
    integer(c_int)               ,intent(inout)      ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! bool Filled() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_Filled ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_Filled ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_Filled')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool StorageOptimized() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_StorageOptimized ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_StorageOptimized ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_StorageOptimized')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool IndicesAreGlobal() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_IndicesAreGlobal ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_IndicesAreGlobal ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_IndicesAreGlobal')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool IndicesAreLocal() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_IndicesAreLocal ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_IndicesAreLocal ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_IndicesAreLocal')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool LowerTriangular() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_LowerTriangular ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_LowerTriangular ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_LowerTriangular')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool UpperTriangular() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_UpperTriangular ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_UpperTriangular ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_UpperTriangular')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool NoDiagonal() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_NoDiagonal ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NoDiagonal ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NoDiagonal')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool MyGlobalRow(int GID) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_MyGlobalRow ( CT_Epetra_CrsGraph_ID_t selfID, int GID );

  function Epetra_CrsGraph_MyGlobalRow ( selfID, GID ) result(that) &
        bind(C,name='Epetra_CrsGraph_MyGlobalRow')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t ,c_int
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GID
  end function


  !> <BR> Original C++ prototype:
  !! bool HaveColMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_HaveColMap ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_HaveColMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_HaveColMap')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyRows() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumMyRows ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumMyRows ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumMyRows')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalRows() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumGlobalRows ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumGlobalRows ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumGlobalRows')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyCols() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumMyCols ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumMyCols ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumMyCols')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalCols() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumGlobalCols ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumGlobalCols ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumGlobalCols')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalNonzeros() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumGlobalNonzeros ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumGlobalNonzeros ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumGlobalNonzeros')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalDiagonals() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumGlobalDiagonals ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumGlobalDiagonals ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumGlobalDiagonals')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyDiagonals() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumMyDiagonals ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumMyDiagonals ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumMyDiagonals')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyBlockRows() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumMyBlockRows ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumMyBlockRows ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumMyBlockRows')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalBlockRows() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumGlobalBlockRows ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumGlobalBlockRows ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumGlobalBlockRows')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyBlockCols() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumMyBlockCols ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumMyBlockCols ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumMyBlockCols')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalBlockCols() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumGlobalBlockCols ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumGlobalBlockCols ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumGlobalBlockCols')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyBlockDiagonals() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumMyBlockDiagonals ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumMyBlockDiagonals ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumMyBlockDiagonals')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalBlockDiagonals() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumGlobalBlockDiagonals ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumGlobalBlockDiagonals ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumGlobalBlockDiagonals')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalEntries() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumGlobalEntries ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumGlobalEntries ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumGlobalEntries')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyEntries() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumMyEntries ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumMyEntries ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumMyEntries')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MaxRowDim() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_MaxRowDim ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_MaxRowDim ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_MaxRowDim')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int GlobalMaxRowDim() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_GlobalMaxRowDim ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_GlobalMaxRowDim ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_GlobalMaxRowDim')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MaxColDim() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_MaxColDim ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_MaxColDim ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_MaxColDim')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int GlobalMaxColDim() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_GlobalMaxColDim ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_GlobalMaxColDim ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_GlobalMaxColDim')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyNonzeros() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumMyNonzeros ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_NumMyNonzeros ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumMyNonzeros')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalIndices(int Row) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumGlobalIndices ( CT_Epetra_CrsGraph_ID_t selfID, int Row );

  function Epetra_CrsGraph_NumGlobalIndices ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumGlobalIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int NumAllocatedGlobalIndices(int Row) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumAllocatedGlobalIndices ( CT_Epetra_CrsGraph_ID_t selfID, int Row );

  function Epetra_CrsGraph_NumAllocatedGlobalIndices ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumAllocatedGlobalIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int MaxNumIndices() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_MaxNumIndices ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_MaxNumIndices ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_MaxNumIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int GlobalMaxNumIndices() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_GlobalMaxNumIndices ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_GlobalMaxNumIndices ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_GlobalMaxNumIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MaxNumNonzeros() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_MaxNumNonzeros ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_MaxNumNonzeros ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_MaxNumNonzeros')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int GlobalMaxNumNonzeros() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_GlobalMaxNumNonzeros ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_GlobalMaxNumNonzeros ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_GlobalMaxNumNonzeros')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyIndices(int Row) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumMyIndices ( CT_Epetra_CrsGraph_ID_t selfID, int Row );

  function Epetra_CrsGraph_NumMyIndices ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumMyIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int NumAllocatedMyIndices(int Row) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_NumAllocatedMyIndices ( CT_Epetra_CrsGraph_ID_t selfID, int Row );

  function Epetra_CrsGraph_NumAllocatedMyIndices ( selfID, Row ) result(that) &
        bind(C,name='Epetra_CrsGraph_NumAllocatedMyIndices')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: Row
  end function


  !> <BR> Original C++ prototype:
  !! int IndexBase() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_IndexBase ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_IndexBase ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_IndexBase')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap& RowMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RowMap ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_RowMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_RowMap')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceRowMap(const Epetra_BlockMap& newmap);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_ReplaceRowMap ( CT_Epetra_CrsGraph_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t newmapID );

  function Epetra_CrsGraph_ReplaceRowMap ( selfID, newmapID ) result(that) &
        bind(C,name='Epetra_CrsGraph_ReplaceRowMap')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: newmapID
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceColMap(const Epetra_BlockMap& newmap);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_ReplaceColMap ( CT_Epetra_CrsGraph_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t newmapID );

  function Epetra_CrsGraph_ReplaceColMap ( selfID, newmapID ) result(that) &
        bind(C,name='Epetra_CrsGraph_ReplaceColMap')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: newmapID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap& ColMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ColMap ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_ColMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_ColMap')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap& DomainMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_DomainMap ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_DomainMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_DomainMap')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap& RangeMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_RangeMap ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_RangeMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_RangeMap')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Import* Importer() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Import_ID_t Epetra_CrsGraph_Importer ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_Importer ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_Importer')
    import :: FT_Epetra_Import_ID_t ,FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_Import_ID_t)                                    :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Export* Exporter() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Export_ID_t Epetra_CrsGraph_Exporter ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_Exporter ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_Exporter')
    import :: FT_Epetra_Export_ID_t ,FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_Export_ID_t)                                    :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Comm& Comm() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Epetra_CrsGraph_Comm ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_Comm ( selfID ) result(that) bind(C,name='Epetra_CrsGraph_Comm')
    import :: FT_Epetra_Comm_ID_t ,FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                      :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int LRID(int GRID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_LRID ( CT_Epetra_CrsGraph_ID_t selfID, int GRID_in );

  function Epetra_CrsGraph_LRID ( selfID, GRID_in ) result(that) &
        bind(C,name='Epetra_CrsGraph_LRID')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GRID_in
  end function


  !> <BR> Original C++ prototype:
  !! int GRID(int LRID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_GRID ( CT_Epetra_CrsGraph_ID_t selfID, int LRID_in );

  function Epetra_CrsGraph_GRID ( selfID, LRID_in ) result(that) &
        bind(C,name='Epetra_CrsGraph_GRID')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LRID_in
  end function


  !> <BR> Original C++ prototype:
  !! int LCID(int GCID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_LCID ( CT_Epetra_CrsGraph_ID_t selfID, int GCID_in );

  function Epetra_CrsGraph_LCID ( selfID, GCID_in ) result(that) &
        bind(C,name='Epetra_CrsGraph_LCID')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GCID_in
  end function


  !> <BR> Original C++ prototype:
  !! int GCID(int LCID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_GCID ( CT_Epetra_CrsGraph_ID_t selfID, int LCID_in );

  function Epetra_CrsGraph_GCID ( selfID, LCID_in ) result(that) &
        bind(C,name='Epetra_CrsGraph_GCID')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LCID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyGRID(int GRID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_MyGRID ( CT_Epetra_CrsGraph_ID_t selfID, int GRID_in );

  function Epetra_CrsGraph_MyGRID ( selfID, GRID_in ) result(that) &
        bind(C,name='Epetra_CrsGraph_MyGRID')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t ,c_int
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GRID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyLRID(int LRID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_MyLRID ( CT_Epetra_CrsGraph_ID_t selfID, int LRID_in );

  function Epetra_CrsGraph_MyLRID ( selfID, LRID_in ) result(that) &
        bind(C,name='Epetra_CrsGraph_MyLRID')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t ,c_int
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LRID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyGCID(int GCID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_MyGCID ( CT_Epetra_CrsGraph_ID_t selfID, int GCID_in );

  function Epetra_CrsGraph_MyGCID ( selfID, GCID_in ) result(that) &
        bind(C,name='Epetra_CrsGraph_MyGCID')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t ,c_int
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GCID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyLCID(int LCID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_CrsGraph_MyLCID ( CT_Epetra_CrsGraph_ID_t selfID, int LCID_in );

  function Epetra_CrsGraph_MyLCID ( selfID, LCID_in ) result(that) &
        bind(C,name='Epetra_CrsGraph_MyLCID')
    import :: FT_boolean_t ,FT_Epetra_CrsGraph_ID_t ,c_int
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LCID_in
  end function


  !> <BR> Original C++ prototype:
  !! inline int* operator[]( int Loc ) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_CrsGraph_getRow ( CT_Epetra_CrsGraph_ID_t selfID, int Loc );

  function Epetra_CrsGraph_getRow ( selfID, Loc ) result(that) &
        bind(C,name='Epetra_CrsGraph_getRow')
    import :: c_ptr ,FT_Epetra_CrsGraph_ID_t ,c_int
    
    type(c_ptr)                                                    :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: Loc
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_CrsGraph& operator = (const Epetra_CrsGraph& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_CrsGraph_Assign ( CT_Epetra_CrsGraph_ID_t selfID, 
  !!     CT_Epetra_CrsGraph_ID_t SourceID );

  subroutine Epetra_CrsGraph_Assign ( selfID, SourceID ) &
        bind(C,name='Epetra_CrsGraph_Assign')
    import :: FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: SourceID
  end subroutine


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap& ImportMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_CrsGraph_ImportMap ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_ImportMap ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_ImportMap')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_CrsGraph_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int TransformToLocal();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_TransformToLocal ( CT_Epetra_CrsGraph_ID_t selfID );

  function Epetra_CrsGraph_TransformToLocal ( selfID ) result(that) &
        bind(C,name='Epetra_CrsGraph_TransformToLocal')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int TransformToLocal(const Epetra_BlockMap* DomainMap, const Epetra_BlockMap* RangeMap);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_CrsGraph_TransformToLocal_UsingMaps ( CT_Epetra_CrsGraph_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t DomainMapID, CT_Epetra_BlockMap_ID_t RangeMapID );

  function Epetra_CrsGraph_TransformToLocal_UsingMaps ( selfID, DomainMapID, RangeMapID ) result(that) &
        bind(C,name='Epetra_CrsGraph_TransformToLocal_UsingMaps')
    import :: c_int ,FT_Epetra_CrsGraph_ID_t ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_CrsGraph_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: DomainMapID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: RangeMapID
  end function


!> @}


!> @name Epetra_DistObject interface
!! @{

  ! _________________ Epetra_DistObject interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_DistObject_ID_t Epetra_DistObject_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_DistObject_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_DistObject_Degeneralize')
    import :: FT_Epetra_DistObject_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_DistObject_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_DistObject_Generalize ( CT_Epetra_DistObject_ID_t id );

  function Epetra_DistObject_Generalize ( id ) result(that) &
        bind(C,name='Epetra_DistObject_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_DistObject_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_DistObject_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_DistObject();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_DistObject_Destroy ( CT_Epetra_DistObject_ID_t * selfID );

  subroutine Epetra_DistObject_Destroy ( selfID ) bind(C,name='Epetra_DistObject_Destroy')
    import :: FT_Epetra_DistObject_ID_t
    
    type(FT_Epetra_DistObject_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int Import(const Epetra_SrcDistObject& A, const Epetra_Import& Importer, 
  !!     Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_DistObject_Import ( CT_Epetra_DistObject_ID_t selfID, 
  !!     CT_Epetra_SrcDistObject_ID_t AID, CT_Epetra_Import_ID_t ImporterID, 
  !!     CT_Epetra_CombineMode_E_t CombineMode, CT_Epetra_OffsetIndex_ID_t IndexorID );

  function Epetra_DistObject_Import ( selfID, AID, ImporterID, CombineMode, IndexorID ) result(that) &
        bind(C,name='Epetra_DistObject_Import')
    import :: c_int ,FT_Epetra_DistObject_ID_t ,FT_Epetra_SrcDistObject_ID_t , &
          FT_Epetra_Import_ID_t ,FT_Epetra_CombineMode_E_t ,FT_Epetra_OffsetIndex_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_DistObject_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SrcDistObject_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_Import_ID_t)    ,intent(in)   ,value              :: ImporterID
    integer(FT_Epetra_CombineMode_E_t),intent(in)   ,value              :: CombineMode
    type(FT_Epetra_OffsetIndex_ID_t),intent(in)   ,value              :: IndexorID
  end function


  !> <BR> Original C++ prototype:
  !! int Import(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, 
  !!     Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_DistObject_Import_UsingExporter ( CT_Epetra_DistObject_ID_t selfID, 
  !!     CT_Epetra_SrcDistObject_ID_t AID, CT_Epetra_Export_ID_t ExporterID, 
  !!     CT_Epetra_CombineMode_E_t CombineMode, CT_Epetra_OffsetIndex_ID_t IndexorID );

  function Epetra_DistObject_Import_UsingExporter ( selfID, AID, ExporterID, CombineMode, &
        IndexorID ) result(that) bind(C,name='Epetra_DistObject_Import_UsingExporter')
    import :: c_int ,FT_Epetra_DistObject_ID_t ,FT_Epetra_SrcDistObject_ID_t , &
          FT_Epetra_Export_ID_t ,FT_Epetra_CombineMode_E_t ,FT_Epetra_OffsetIndex_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_DistObject_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SrcDistObject_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_Export_ID_t)    ,intent(in)   ,value              :: ExporterID
    integer(FT_Epetra_CombineMode_E_t),intent(in)   ,value              :: CombineMode
    type(FT_Epetra_OffsetIndex_ID_t),intent(in)   ,value              :: IndexorID
  end function


  !> <BR> Original C++ prototype:
  !! int Export(const Epetra_SrcDistObject& A, const Epetra_Import & Importer, 
  !!     Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_DistObject_Export_UsingImporter ( CT_Epetra_DistObject_ID_t selfID, 
  !!     CT_Epetra_SrcDistObject_ID_t AID, CT_Epetra_Import_ID_t ImporterID, 
  !!     CT_Epetra_CombineMode_E_t CombineMode, CT_Epetra_OffsetIndex_ID_t IndexorID );

  function Epetra_DistObject_Export_UsingImporter ( selfID, AID, ImporterID, CombineMode, &
        IndexorID ) result(that) bind(C,name='Epetra_DistObject_Export_UsingImporter')
    import :: c_int ,FT_Epetra_DistObject_ID_t ,FT_Epetra_SrcDistObject_ID_t , &
          FT_Epetra_Import_ID_t ,FT_Epetra_CombineMode_E_t ,FT_Epetra_OffsetIndex_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_DistObject_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SrcDistObject_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_Import_ID_t)    ,intent(in)   ,value              :: ImporterID
    integer(FT_Epetra_CombineMode_E_t),intent(in)   ,value              :: CombineMode
    type(FT_Epetra_OffsetIndex_ID_t),intent(in)   ,value              :: IndexorID
  end function


  !> <BR> Original C++ prototype:
  !! int Export(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, 
  !!     Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_DistObject_Export ( CT_Epetra_DistObject_ID_t selfID, 
  !!     CT_Epetra_SrcDistObject_ID_t AID, CT_Epetra_Export_ID_t ExporterID, 
  !!     CT_Epetra_CombineMode_E_t CombineMode, CT_Epetra_OffsetIndex_ID_t IndexorID );

  function Epetra_DistObject_Export ( selfID, AID, ExporterID, CombineMode, IndexorID ) result(that) &
        bind(C,name='Epetra_DistObject_Export')
    import :: c_int ,FT_Epetra_DistObject_ID_t ,FT_Epetra_SrcDistObject_ID_t , &
          FT_Epetra_Export_ID_t ,FT_Epetra_CombineMode_E_t ,FT_Epetra_OffsetIndex_ID_t
    
    integer(c_int)                                                   :: that
    type(FT_Epetra_DistObject_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SrcDistObject_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_Export_ID_t)    ,intent(in)   ,value              :: ExporterID
    integer(FT_Epetra_CombineMode_E_t),intent(in)   ,value              :: CombineMode
    type(FT_Epetra_OffsetIndex_ID_t),intent(in)   ,value              :: IndexorID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap& Map() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_DistObject_Map ( CT_Epetra_DistObject_ID_t selfID );

  function Epetra_DistObject_Map ( selfID ) result(that) &
        bind(C,name='Epetra_DistObject_Map')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_DistObject_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                    :: that
    type(FT_Epetra_DistObject_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Comm& Comm() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Epetra_DistObject_Comm ( CT_Epetra_DistObject_ID_t selfID );

  function Epetra_DistObject_Comm ( selfID ) result(that) &
        bind(C,name='Epetra_DistObject_Comm')
    import :: FT_Epetra_Comm_ID_t ,FT_Epetra_DistObject_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                        :: that
    type(FT_Epetra_DistObject_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool DistributedGlobal() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_DistObject_DistributedGlobal ( CT_Epetra_DistObject_ID_t selfID );

  function Epetra_DistObject_DistributedGlobal ( selfID ) result(that) &
        bind(C,name='Epetra_DistObject_DistributedGlobal')
    import :: FT_boolean_t ,FT_Epetra_DistObject_ID_t
    
    integer(FT_boolean_t)                                            :: that
    type(FT_Epetra_DistObject_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_Vector interface
!! @{

  ! _________________ Epetra_Vector interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Vector_ID_t Epetra_Vector_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Vector_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Vector_Degeneralize')
    import :: FT_Epetra_Vector_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Vector_ID_t)                                   :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Vector_Generalize ( CT_Epetra_Vector_ID_t id );

  function Epetra_Vector_Generalize ( id ) result(that) &
        bind(C,name='Epetra_Vector_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Vector_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Vector(const Epetra_BlockMap& Map, bool zeroOut = true);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Vector_ID_t Epetra_Vector_Create ( CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut );

  function Epetra_Vector_Create ( MapID, zeroOut ) result(that) &
        bind(C,name='Epetra_Vector_Create')
    import :: FT_Epetra_Vector_ID_t ,FT_Epetra_BlockMap_ID_t ,FT_boolean_t
    
    type(FT_Epetra_Vector_ID_t)                                   :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: MapID
    integer(FT_boolean_t)       ,intent(in)   ,value              :: zeroOut
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Vector(const Epetra_Vector& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( CT_Epetra_Vector_ID_t SourceID );

  function Epetra_Vector_Duplicate ( SourceID ) result(that) &
        bind(C,name='Epetra_Vector_Duplicate')
    import :: FT_Epetra_Vector_ID_t
    
    type(FT_Epetra_Vector_ID_t)                                   :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: SourceID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Vector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double *V);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Vector_ID_t Epetra_Vector_Create_FromArray ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_BlockMap_ID_t MapID, double * V );

  function Epetra_Vector_Create_FromArray ( CV, MapID, V ) result(that) &
        bind(C,name='Epetra_Vector_Create_FromArray')
    import :: FT_Epetra_Vector_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_BlockMap_ID_t , &
          c_double
    
    type(FT_Epetra_Vector_ID_t)                                   :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: MapID
    real(c_double)                                  ,dimension(*) :: V
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Vector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int Index);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Vector_ID_t Epetra_Vector_FromSource ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_MultiVector_ID_t SourceID, int Index );

  function Epetra_Vector_FromSource ( CV, SourceID, Index ) result(that) &
        bind(C,name='Epetra_Vector_FromSource')
    import :: FT_Epetra_Vector_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_MultiVector_ID_t , &
          c_int
    
    type(FT_Epetra_Vector_ID_t)                                   :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: SourceID
    integer(c_int)              ,intent(in)   ,value              :: Index
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Vector ();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Vector_Destroy ( CT_Epetra_Vector_ID_t * selfID );

  subroutine Epetra_Vector_Destroy ( selfID ) bind(C,name='Epetra_Vector_Destroy')
    import :: FT_Epetra_Vector_ID_t
    
    type(FT_Epetra_Vector_ID_t)                                   :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValues(int NumEntries, double * Values, int * Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_ReplaceGlobalValues ( CT_Epetra_Vector_ID_t selfID, int NumEntries, 
  !!     double * Values, int * Indices );

  function Epetra_Vector_ReplaceGlobalValues ( selfID, NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_Vector_ReplaceGlobalValues')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumEntries
    real(c_double)                                  ,dimension(*) :: Values
    integer(c_int)                                  ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceMyValues(int NumEntries, double * Values, int * Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_ReplaceMyValues ( CT_Epetra_Vector_ID_t selfID, int NumEntries, 
  !!     double * Values, int * Indices );

  function Epetra_Vector_ReplaceMyValues ( selfID, NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_Vector_ReplaceMyValues')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumEntries
    real(c_double)                                  ,dimension(*) :: Values
    integer(c_int)                                  ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValues(int NumEntries, double * Values, int * Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_SumIntoGlobalValues ( CT_Epetra_Vector_ID_t selfID, int NumEntries, 
  !!     double * Values, int * Indices );

  function Epetra_Vector_SumIntoGlobalValues ( selfID, NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_Vector_SumIntoGlobalValues')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumEntries
    real(c_double)                                  ,dimension(*) :: Values
    integer(c_int)                                  ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoMyValues(int NumEntries, double * Values, int * Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_SumIntoMyValues ( CT_Epetra_Vector_ID_t selfID, int NumEntries, 
  !!     double * Values, int * Indices );

  function Epetra_Vector_SumIntoMyValues ( selfID, NumEntries, Values, Indices ) result(that) &
        bind(C,name='Epetra_Vector_SumIntoMyValues')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumEntries
    real(c_double)                                  ,dimension(*) :: Values
    integer(c_int)                                  ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_ReplaceGlobalValues_BlockPos ( CT_Epetra_Vector_ID_t selfID, int NumEntries, 
  !!     int BlockOffset, double * Values, int * Indices );

  function Epetra_Vector_ReplaceGlobalValues_BlockPos ( selfID, NumEntries, BlockOffset, &
        Values, Indices ) result(that) &
        bind(C,name='Epetra_Vector_ReplaceGlobalValues_BlockPos')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumEntries
    integer(c_int)              ,intent(in)   ,value              :: BlockOffset
    real(c_double)                                  ,dimension(*) :: Values
    integer(c_int)                                  ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_ReplaceMyValues_BlockPos ( CT_Epetra_Vector_ID_t selfID, int NumEntries, 
  !!     int BlockOffset, double * Values, int * Indices );

  function Epetra_Vector_ReplaceMyValues_BlockPos ( selfID, NumEntries, BlockOffset, Values, &
        Indices ) result(that) bind(C,name='Epetra_Vector_ReplaceMyValues_BlockPos')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumEntries
    integer(c_int)              ,intent(in)   ,value              :: BlockOffset
    real(c_double)                                  ,dimension(*) :: Values
    integer(c_int)                                  ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_SumIntoGlobalValues_BlockPos ( CT_Epetra_Vector_ID_t selfID, int NumEntries, 
  !!     int BlockOffset, double * Values, int * Indices );

  function Epetra_Vector_SumIntoGlobalValues_BlockPos ( selfID, NumEntries, BlockOffset, &
        Values, Indices ) result(that) &
        bind(C,name='Epetra_Vector_SumIntoGlobalValues_BlockPos')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumEntries
    integer(c_int)              ,intent(in)   ,value              :: BlockOffset
    real(c_double)                                  ,dimension(*) :: Values
    integer(c_int)                                  ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_SumIntoMyValues_BlockPos ( CT_Epetra_Vector_ID_t selfID, int NumEntries, 
  !!     int BlockOffset, double * Values, int * Indices );

  function Epetra_Vector_SumIntoMyValues_BlockPos ( selfID, NumEntries, BlockOffset, Values, &
        Indices ) result(that) bind(C,name='Epetra_Vector_SumIntoMyValues_BlockPos')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: NumEntries
    integer(c_int)              ,intent(in)   ,value              :: BlockOffset
    real(c_double)                                  ,dimension(*) :: Values
    integer(c_int)                                  ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractCopy(double *V) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_ExtractCopy ( CT_Epetra_Vector_ID_t selfID, double * V );

  function Epetra_Vector_ExtractCopy ( selfID, V ) result(that) &
        bind(C,name='Epetra_Vector_ExtractCopy')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: V
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractView(double **V) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Vector_ExtractView ( CT_Epetra_Vector_ID_t selfID, double ** V );

  function Epetra_Vector_ExtractView ( selfID, V ) result(that) &
        bind(C,name='Epetra_Vector_ExtractView')
    import :: c_int ,FT_Epetra_Vector_ID_t ,c_double
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    real(c_double)                                  ,dimension(*) :: V
  end function


  !> <BR> Original C++ prototype:
  !! const double& operator [] (int index) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_Vector_getElement ( CT_Epetra_Vector_ID_t selfID, int index );

  function Epetra_Vector_getElement ( selfID, index ) result(that) &
        bind(C,name='Epetra_Vector_getElement')
    import :: c_double ,FT_Epetra_Vector_ID_t ,c_int
    
    real(c_double)                                                :: that
    type(FT_Epetra_Vector_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: index
  end function


!> @}


!> @name Epetra_Export interface
!! @{

  ! _________________ Epetra_Export interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Export_ID_t Epetra_Export_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Export_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Export_Degeneralize')
    import :: FT_Epetra_Export_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Export_ID_t)                                   :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Export_Generalize ( CT_Epetra_Export_ID_t id );

  function Epetra_Export_Generalize ( id ) result(that) &
        bind(C,name='Epetra_Export_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Export_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Export( const Epetra_BlockMap & SourceMap, const Epetra_BlockMap & TargetMap );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Export_ID_t Epetra_Export_Create ( CT_Epetra_BlockMap_ID_t SourceMapID, 
  !!     CT_Epetra_BlockMap_ID_t TargetMapID );

  function Epetra_Export_Create ( SourceMapID, TargetMapID ) result(that) &
        bind(C,name='Epetra_Export_Create')
    import :: FT_Epetra_Export_ID_t ,FT_Epetra_BlockMap_ID_t
    
    type(FT_Epetra_Export_ID_t)                                   :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: SourceMapID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: TargetMapID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Export(const Epetra_Export& Exporter);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( CT_Epetra_Export_ID_t ExporterID );

  function Epetra_Export_Duplicate ( ExporterID ) result(that) &
        bind(C,name='Epetra_Export_Duplicate')
    import :: FT_Epetra_Export_ID_t
    
    type(FT_Epetra_Export_ID_t)                                   :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: ExporterID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Export(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID );

  subroutine Epetra_Export_Destroy ( selfID ) bind(C,name='Epetra_Export_Destroy')
    import :: FT_Epetra_Export_ID_t
    
    type(FT_Epetra_Export_ID_t)                                   :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int NumSameIDs() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Export_NumSameIDs ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_NumSameIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Export_NumSameIDs')
    import :: c_int ,FT_Epetra_Export_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumPermuteIDs() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Export_NumPermuteIDs ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_NumPermuteIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Export_NumPermuteIDs')
    import :: c_int ,FT_Epetra_Export_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * PermuteFromLIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Export_PermuteFromLIDs ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_PermuteFromLIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Export_PermuteFromLIDs')
    import :: c_ptr ,FT_Epetra_Export_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * PermuteToLIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Export_PermuteToLIDs ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_PermuteToLIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Export_PermuteToLIDs')
    import :: c_ptr ,FT_Epetra_Export_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumRemoteIDs() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Export_NumRemoteIDs ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_NumRemoteIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Export_NumRemoteIDs')
    import :: c_int ,FT_Epetra_Export_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * RemoteLIDs() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Export_RemoteLIDs ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_RemoteLIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Export_RemoteLIDs')
    import :: c_ptr ,FT_Epetra_Export_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumExportIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Export_NumExportIDs ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_NumExportIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Export_NumExportIDs')
    import :: c_int ,FT_Epetra_Export_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * ExportLIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Export_ExportLIDs ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_ExportLIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Export_ExportLIDs')
    import :: c_ptr ,FT_Epetra_Export_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * ExportPIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Export_ExportPIDs ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_ExportPIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Export_ExportPIDs')
    import :: c_ptr ,FT_Epetra_Export_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumSend() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Export_NumSend ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_NumSend ( selfID ) result(that) &
        bind(C,name='Epetra_Export_NumSend')
    import :: c_int ,FT_Epetra_Export_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumRecv() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Export_NumRecv ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_NumRecv ( selfID ) result(that) &
        bind(C,name='Epetra_Export_NumRecv')
    import :: c_int ,FT_Epetra_Export_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap & SourceMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_Export_SourceMap ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_SourceMap ( selfID ) result(that) &
        bind(C,name='Epetra_Export_SourceMap')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_Export_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap & TargetMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_TargetMap ( selfID ) result(that) &
        bind(C,name='Epetra_Export_TargetMap')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_Export_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Distributor & Distributor() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Distributor_ID_t Epetra_Export_Distributor ( CT_Epetra_Export_ID_t selfID );

  function Epetra_Export_Distributor ( selfID ) result(that) &
        bind(C,name='Epetra_Export_Distributor')
    import :: FT_Epetra_Distributor_ID_t ,FT_Epetra_Export_ID_t
    
    type(FT_Epetra_Distributor_ID_t)                                  :: that
    type(FT_Epetra_Export_ID_t) ,intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_Map interface
!! @{

  ! _________________ Epetra_Map interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_Map_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Map_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Map_Degeneralize')
    import :: FT_Epetra_Map_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Map_ID_t)                                      :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Map_Generalize ( CT_Epetra_Map_ID_t id );

  function Epetra_Map_Generalize ( id ) result(that) bind(C,name='Epetra_Map_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Map_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Map_ID_t)    ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_Map_Create ( int NumGlobalElements, int IndexBase, 
  !!     CT_Epetra_Comm_ID_t CommID );

  function Epetra_Map_Create ( NumGlobalElements, IndexBase, CommID ) result(that) &
        bind(C,name='Epetra_Map_Create')
    import :: FT_Epetra_Map_ID_t ,c_int ,FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_Map_ID_t)                                      :: that
    integer(c_int)              ,intent(in)   ,value              :: NumGlobalElements
    integer(c_int)              ,intent(in)   ,value              :: IndexBase
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( int NumGlobalElements, int NumMyElements, 
  !!     int IndexBase, CT_Epetra_Comm_ID_t CommID );

  function Epetra_Map_Create_Linear ( NumGlobalElements, NumMyElements, IndexBase, CommID ) result(that) &
        bind(C,name='Epetra_Map_Create_Linear')
    import :: FT_Epetra_Map_ID_t ,c_int ,FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_Map_ID_t)                                      :: that
    integer(c_int)              ,intent(in)   ,value              :: NumGlobalElements
    integer(c_int)              ,intent(in)   ,value              :: NumMyElements
    integer(c_int)              ,intent(in)   ,value              :: IndexBase
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Map(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, 
  !!     int IndexBase, const Epetra_Comm& Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( int NumGlobalElements, int NumMyElements, 
  !!     const int * MyGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  function Epetra_Map_Create_Arbitrary ( NumGlobalElements, NumMyElements, MyGlobalElements, &
        IndexBase, CommID ) result(that) bind(C,name='Epetra_Map_Create_Arbitrary')
    import :: FT_Epetra_Map_ID_t ,c_int ,FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_Map_ID_t)                                      :: that
    integer(c_int)              ,intent(in)   ,value              :: NumGlobalElements
    integer(c_int)              ,intent(in)   ,value              :: NumMyElements
    integer(c_int)              ,intent(in)         ,dimension(*) :: MyGlobalElements
    integer(c_int)              ,intent(in)   ,value              :: IndexBase
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Map(const Epetra_Map& map);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );

  function Epetra_Map_Duplicate ( mapID ) result(that) bind(C,name='Epetra_Map_Duplicate')
    import :: FT_Epetra_Map_ID_t
    
    type(FT_Epetra_Map_ID_t)                                      :: that
    type(FT_Epetra_Map_ID_t)    ,intent(in)   ,value              :: mapID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Map(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID );

  subroutine Epetra_Map_Destroy ( selfID ) bind(C,name='Epetra_Map_Destroy')
    import :: FT_Epetra_Map_ID_t
    
    type(FT_Epetra_Map_ID_t)                                      :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! Epetra_Map & operator=(const Epetra_Map & map);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Map_Assign ( CT_Epetra_Map_ID_t selfID, CT_Epetra_Map_ID_t mapID );

  subroutine Epetra_Map_Assign ( selfID, mapID ) bind(C,name='Epetra_Map_Assign')
    import :: FT_Epetra_Map_ID_t
    
    type(FT_Epetra_Map_ID_t)    ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Map_ID_t)    ,intent(in)   ,value              :: mapID
  end subroutine


!> @}


!> @name Epetra_BlockMap interface
!! @{

  ! _________________ Epetra_BlockMap interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_BlockMap_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_BlockMap_Degeneralize')
    import :: FT_Epetra_BlockMap_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_BlockMap_Generalize ( CT_Epetra_BlockMap_ID_t id );

  function Epetra_BlockMap_Generalize ( id ) result(that) &
        bind(C,name='Epetra_BlockMap_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_BlockMap_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, 
  !!     const Epetra_Comm& Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( int NumGlobalElements, int ElementSize, 
  !!     int IndexBase, CT_Epetra_Comm_ID_t CommID );

  function Epetra_BlockMap_Create ( NumGlobalElements, ElementSize, IndexBase, CommID ) result(that) &
        bind(C,name='Epetra_BlockMap_Create')
    import :: FT_Epetra_BlockMap_ID_t ,c_int ,FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    integer(c_int)               ,intent(in)   ,value              :: NumGlobalElements
    integer(c_int)               ,intent(in)   ,value              :: ElementSize
    integer(c_int)               ,intent(in)   ,value              :: IndexBase
    type(FT_Epetra_Comm_ID_t)    ,intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int ElementSize, int IndexBase, 
  !!     const Epetra_Comm& Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( int NumGlobalElements, 
  !!     int NumMyElements, int ElementSize, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  function Epetra_BlockMap_Create_Linear ( NumGlobalElements, NumMyElements, ElementSize, &
        IndexBase, CommID ) result(that) bind(C,name='Epetra_BlockMap_Create_Linear')
    import :: FT_Epetra_BlockMap_ID_t ,c_int ,FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    integer(c_int)               ,intent(in)   ,value              :: NumGlobalElements
    integer(c_int)               ,intent(in)   ,value              :: NumMyElements
    integer(c_int)               ,intent(in)   ,value              :: ElementSize
    integer(c_int)               ,intent(in)   ,value              :: IndexBase
    type(FT_Epetra_Comm_ID_t)    ,intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, 
  !!     int ElementSize, int IndexBase, const Epetra_Comm& Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( int NumGlobalElements, 
  !!     int NumMyElements, const int * MyGlobalElements, int ElementSize, int IndexBase, 
  !!     CT_Epetra_Comm_ID_t CommID );

  function Epetra_BlockMap_Create_Arbitrary ( NumGlobalElements, NumMyElements, &
        MyGlobalElements, ElementSize, IndexBase, CommID ) result(that) &
        bind(C,name='Epetra_BlockMap_Create_Arbitrary')
    import :: FT_Epetra_BlockMap_ID_t ,c_int ,FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    integer(c_int)               ,intent(in)   ,value              :: NumGlobalElements
    integer(c_int)               ,intent(in)   ,value              :: NumMyElements
    integer(c_int)               ,intent(in)         ,dimension(*) :: MyGlobalElements
    integer(c_int)               ,intent(in)   ,value              :: ElementSize
    integer(c_int)               ,intent(in)   ,value              :: IndexBase
    type(FT_Epetra_Comm_ID_t)    ,intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, 
  !!     const int *ElementSizeList, int IndexBase, const Epetra_Comm& Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( int NumGlobalElements, 
  !!     int NumMyElements, const int * MyGlobalElements, const int * ElementSizeList, 
  !!     int IndexBase, CT_Epetra_Comm_ID_t CommID );

  function Epetra_BlockMap_Create_Variable ( NumGlobalElements, NumMyElements, &
        MyGlobalElements, ElementSizeList, IndexBase, CommID ) result(that) &
        bind(C,name='Epetra_BlockMap_Create_Variable')
    import :: FT_Epetra_BlockMap_ID_t ,c_int ,FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    integer(c_int)               ,intent(in)   ,value              :: NumGlobalElements
    integer(c_int)               ,intent(in)   ,value              :: NumMyElements
    integer(c_int)               ,intent(in)         ,dimension(*) :: MyGlobalElements
    integer(c_int)               ,intent(in)         ,dimension(*) :: ElementSizeList
    integer(c_int)               ,intent(in)   ,value              :: IndexBase
    type(FT_Epetra_Comm_ID_t)    ,intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_BlockMap(const Epetra_BlockMap& map);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( CT_Epetra_BlockMap_ID_t mapID );

  function Epetra_BlockMap_Duplicate ( mapID ) result(that) &
        bind(C,name='Epetra_BlockMap_Duplicate')
    import :: FT_Epetra_BlockMap_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: mapID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_BlockMap(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BlockMap_Destroy ( CT_Epetra_BlockMap_ID_t * selfID );

  subroutine Epetra_BlockMap_Destroy ( selfID ) bind(C,name='Epetra_BlockMap_Destroy')
    import :: FT_Epetra_BlockMap_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_RemoteIDList ( CT_Epetra_BlockMap_ID_t selfID, int NumIDs, 
  !!     const int * GIDList, int * PIDList, int * LIDList );

  function Epetra_BlockMap_RemoteIDList ( selfID, NumIDs, GIDList, PIDList, LIDList ) result(that) &
        bind(C,name='Epetra_BlockMap_RemoteIDList')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: NumIDs
    integer(c_int)               ,intent(in)         ,dimension(*) :: GIDList
    integer(c_int)                                   ,dimension(*) :: PIDList
    integer(c_int)                                   ,dimension(*) :: LIDList
  end function


  !> <BR> Original C++ prototype:
  !! int RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList, 
  !!     int * SizeList) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_RemoteIDList_WithSize ( CT_Epetra_BlockMap_ID_t selfID, int NumIDs, 
  !!     const int * GIDList, int * PIDList, int * LIDList, int * SizeList );

  function Epetra_BlockMap_RemoteIDList_WithSize ( selfID, NumIDs, GIDList, PIDList, &
        LIDList, SizeList ) result(that) &
        bind(C,name='Epetra_BlockMap_RemoteIDList_WithSize')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: NumIDs
    integer(c_int)               ,intent(in)         ,dimension(*) :: GIDList
    integer(c_int)                                   ,dimension(*) :: PIDList
    integer(c_int)                                   ,dimension(*) :: LIDList
    integer(c_int)                                   ,dimension(*) :: SizeList
  end function


  !> <BR> Original C++ prototype:
  !! int LID(int GID) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_LID ( CT_Epetra_BlockMap_ID_t selfID, int GID );

  function Epetra_BlockMap_LID ( selfID, GID ) result(that) &
        bind(C,name='Epetra_BlockMap_LID')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GID
  end function


  !> <BR> Original C++ prototype:
  !! int GID(int LID) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_GID ( CT_Epetra_BlockMap_ID_t selfID, int LID );

  function Epetra_BlockMap_GID ( selfID, LID ) result(that) &
        bind(C,name='Epetra_BlockMap_GID')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LID
  end function


  !> <BR> Original C++ prototype:
  !! int FindLocalElementID(int PointID, int & ElementID, int & ElementOffset) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_FindLocalElementID ( CT_Epetra_BlockMap_ID_t selfID, int PointID, 
  !!     int * ElementID, int * ElementOffset );

  function Epetra_BlockMap_FindLocalElementID ( selfID, PointID, ElementID, ElementOffset ) result(that) &
        bind(C,name='Epetra_BlockMap_FindLocalElementID')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: PointID
    integer(c_int)               ,intent(inout)                    :: ElementID
    integer(c_int)               ,intent(inout)                    :: ElementOffset
  end function


  !> <BR> Original C++ prototype:
  !! bool MyGID(int GID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_BlockMap_MyGID ( CT_Epetra_BlockMap_ID_t selfID, int GID_in );

  function Epetra_BlockMap_MyGID ( selfID, GID_in ) result(that) &
        bind(C,name='Epetra_BlockMap_MyGID')
    import :: FT_boolean_t ,FT_Epetra_BlockMap_ID_t ,c_int
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: GID_in
  end function


  !> <BR> Original C++ prototype:
  !! bool MyLID(int LID_in) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_BlockMap_MyLID ( CT_Epetra_BlockMap_ID_t selfID, int LID_in );

  function Epetra_BlockMap_MyLID ( selfID, LID_in ) result(that) &
        bind(C,name='Epetra_BlockMap_MyLID')
    import :: FT_boolean_t ,FT_Epetra_BlockMap_ID_t ,c_int
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LID_in
  end function


  !> <BR> Original C++ prototype:
  !! int MinAllGID() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MinAllGID ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MinAllGID ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MinAllGID')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MaxAllGID() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MaxAllGID ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MaxAllGID ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MaxAllGID')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MinMyGID() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MinMyGID ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MinMyGID ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MinMyGID')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MaxMyGID() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MaxMyGID ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MaxMyGID ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MaxMyGID')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MinLID() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MinLID ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MinLID ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MinLID')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MaxLID() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MaxLID ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MaxLID ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MaxLID')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalElements() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_NumGlobalElements ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_NumGlobalElements ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_NumGlobalElements')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyElements() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_NumMyElements ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_NumMyElements ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_NumMyElements')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MyGlobalElements(int * MyGlobalElementList) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MyGlobalElements_Fill ( CT_Epetra_BlockMap_ID_t selfID, 
  !!     int * MyGlobalElementList );

  function Epetra_BlockMap_MyGlobalElements_Fill ( selfID, MyGlobalElementList ) result(that) &
        bind(C,name='Epetra_BlockMap_MyGlobalElements_Fill')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                   ,dimension(*) :: MyGlobalElementList
  end function


  !> <BR> Original C++ prototype:
  !! int ElementSize() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_ElementSize_Const ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_ElementSize_Const ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_ElementSize_Const')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ElementSize(int LID) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_ElementSize ( CT_Epetra_BlockMap_ID_t selfID, int LID );

  function Epetra_BlockMap_ElementSize ( selfID, LID ) result(that) &
        bind(C,name='Epetra_BlockMap_ElementSize')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LID
  end function


  !> <BR> Original C++ prototype:
  !! int FirstPointInElement(int LID) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_FirstPointInElement ( CT_Epetra_BlockMap_ID_t selfID, int LID );

  function Epetra_BlockMap_FirstPointInElement ( selfID, LID ) result(that) &
        bind(C,name='Epetra_BlockMap_FirstPointInElement')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)               ,intent(in)   ,value              :: LID
  end function


  !> <BR> Original C++ prototype:
  !! int IndexBase() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_IndexBase ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_IndexBase ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_IndexBase')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumGlobalPoints() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_NumGlobalPoints ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_NumGlobalPoints ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_NumGlobalPoints')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyPoints() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_NumMyPoints ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_NumMyPoints ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_NumMyPoints')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MinMyElementSize() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MinMyElementSize ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MinMyElementSize ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MinMyElementSize')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MaxMyElementSize() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MaxMyElementSize ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MaxMyElementSize ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MaxMyElementSize')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MinElementSize() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MinElementSize ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MinElementSize ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MinElementSize')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int MaxElementSize() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_MaxElementSize ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MaxElementSize ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MaxElementSize')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool UniqueGIDs() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_BlockMap_UniqueGIDs ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_UniqueGIDs ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_UniqueGIDs')
    import :: FT_boolean_t ,FT_Epetra_BlockMap_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool ConstantElementSize() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_BlockMap_ConstantElementSize ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_ConstantElementSize ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_ConstantElementSize')
    import :: FT_boolean_t ,FT_Epetra_BlockMap_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool SameAs(const Epetra_BlockMap & Map) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_BlockMap_SameAs ( CT_Epetra_BlockMap_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t MapID );

  function Epetra_BlockMap_SameAs ( selfID, MapID ) result(that) &
        bind(C,name='Epetra_BlockMap_SameAs')
    import :: FT_boolean_t ,FT_Epetra_BlockMap_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: MapID
  end function


  !> <BR> Original C++ prototype:
  !! bool PointSameAs(const Epetra_BlockMap & Map) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_BlockMap_PointSameAs ( CT_Epetra_BlockMap_ID_t selfID, 
  !!     CT_Epetra_BlockMap_ID_t MapID );

  function Epetra_BlockMap_PointSameAs ( selfID, MapID ) result(that) &
        bind(C,name='Epetra_BlockMap_PointSameAs')
    import :: FT_boolean_t ,FT_Epetra_BlockMap_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: MapID
  end function


  !> <BR> Original C++ prototype:
  !! bool LinearMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_BlockMap_LinearMap ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_LinearMap ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_LinearMap')
    import :: FT_boolean_t ,FT_Epetra_BlockMap_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool DistributedGlobal() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_BlockMap_DistributedGlobal ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_DistributedGlobal ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_DistributedGlobal')
    import :: FT_boolean_t ,FT_Epetra_BlockMap_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * MyGlobalElements() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_BlockMap_MyGlobalElements ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_MyGlobalElements ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_MyGlobalElements')
    import :: c_ptr ,FT_Epetra_BlockMap_ID_t
    
    type(c_ptr)                                                    :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * FirstPointInElementList() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_BlockMap_FirstPointInElementList ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_FirstPointInElementList ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_FirstPointInElementList')
    import :: c_ptr ,FT_Epetra_BlockMap_ID_t
    
    type(c_ptr)                                                    :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * ElementSizeList() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_BlockMap_ElementSizeList ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_ElementSizeList ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_ElementSizeList')
    import :: c_ptr ,FT_Epetra_BlockMap_ID_t
    
    type(c_ptr)                                                    :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * PointToElementList() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_BlockMap_PointToElementList ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_PointToElementList ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_PointToElementList')
    import :: c_ptr ,FT_Epetra_BlockMap_ID_t
    
    type(c_ptr)                                                    :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int ElementSizeList(int * ElementSizeList)const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_ElementSizeList_Fill ( CT_Epetra_BlockMap_ID_t selfID, 
  !!     int * ElementSizeList );

  function Epetra_BlockMap_ElementSizeList_Fill ( selfID, ElementSizeList ) result(that) &
        bind(C,name='Epetra_BlockMap_ElementSizeList_Fill')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                   ,dimension(*) :: ElementSizeList
  end function


  !> <BR> Original C++ prototype:
  !! int FirstPointInElementList(int * FirstPointInElementList)const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_FirstPointInElementList_Fill ( CT_Epetra_BlockMap_ID_t selfID, 
  !!     int * FirstPointInElementList );

  function Epetra_BlockMap_FirstPointInElementList_Fill ( selfID, FirstPointInElementList ) result(that) &
        bind(C,name='Epetra_BlockMap_FirstPointInElementList_Fill')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                   ,dimension(*) :: FirstPointInElementList
  end function


  !> <BR> Original C++ prototype:
  !! int PointToElementList(int * PointToElementList) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_BlockMap_PointToElementList_Fill ( CT_Epetra_BlockMap_ID_t selfID, 
  !!     int * PointToElementList );

  function Epetra_BlockMap_PointToElementList_Fill ( selfID, PointToElementList ) result(that) &
        bind(C,name='Epetra_BlockMap_PointToElementList_Fill')
    import :: c_int ,FT_Epetra_BlockMap_ID_t
    
    integer(c_int)                                                 :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                   ,dimension(*) :: PointToElementList
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_Comm & Comm() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Comm_ID_t Epetra_BlockMap_Comm ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_Comm ( selfID ) result(that) bind(C,name='Epetra_BlockMap_Comm')
    import :: FT_Epetra_Comm_ID_t ,FT_Epetra_BlockMap_ID_t
    
    type(FT_Epetra_Comm_ID_t)                                      :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool IsOneToOne() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_BlockMap_IsOneToOne ( CT_Epetra_BlockMap_ID_t selfID );

  function Epetra_BlockMap_IsOneToOne ( selfID ) result(that) &
        bind(C,name='Epetra_BlockMap_IsOneToOne')
    import :: FT_boolean_t ,FT_Epetra_BlockMap_ID_t
    
    integer(FT_boolean_t)                                          :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_BlockMap & operator=(const Epetra_BlockMap & map);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_BlockMap_Assign ( CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID );

  subroutine Epetra_BlockMap_Assign ( selfID, mapID ) bind(C,name='Epetra_BlockMap_Assign')
    import :: FT_Epetra_BlockMap_ID_t
    
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: mapID
  end subroutine


!> @}


!> @name Epetra_Import interface
!! @{

  ! _________________ Epetra_Import interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Import_ID_t Epetra_Import_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Import_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Import_Degeneralize')
    import :: FT_Epetra_Import_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Import_ID_t)                                   :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Import_Generalize ( CT_Epetra_Import_ID_t id );

  function Epetra_Import_Generalize ( id ) result(that) &
        bind(C,name='Epetra_Import_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Import_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Import( const Epetra_BlockMap & TargetMap, const Epetra_BlockMap & SourceMap );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Import_ID_t Epetra_Import_Create ( CT_Epetra_BlockMap_ID_t TargetMapID, 
  !!     CT_Epetra_BlockMap_ID_t SourceMapID );

  function Epetra_Import_Create ( TargetMapID, SourceMapID ) result(that) &
        bind(C,name='Epetra_Import_Create')
    import :: FT_Epetra_Import_ID_t ,FT_Epetra_BlockMap_ID_t
    
    type(FT_Epetra_Import_ID_t)                                   :: that
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: TargetMapID
    type(FT_Epetra_BlockMap_ID_t),intent(in)   ,value              :: SourceMapID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Import(const Epetra_Import& Importer);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( CT_Epetra_Import_ID_t ImporterID );

  function Epetra_Import_Duplicate ( ImporterID ) result(that) &
        bind(C,name='Epetra_Import_Duplicate')
    import :: FT_Epetra_Import_ID_t
    
    type(FT_Epetra_Import_ID_t)                                   :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: ImporterID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Import(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Import_Destroy ( CT_Epetra_Import_ID_t * selfID );

  subroutine Epetra_Import_Destroy ( selfID ) bind(C,name='Epetra_Import_Destroy')
    import :: FT_Epetra_Import_ID_t
    
    type(FT_Epetra_Import_ID_t)                                   :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int NumSameIDs() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Import_NumSameIDs ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_NumSameIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Import_NumSameIDs')
    import :: c_int ,FT_Epetra_Import_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumPermuteIDs() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Import_NumPermuteIDs ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_NumPermuteIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Import_NumPermuteIDs')
    import :: c_int ,FT_Epetra_Import_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * PermuteFromLIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Import_PermuteFromLIDs ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_PermuteFromLIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Import_PermuteFromLIDs')
    import :: c_ptr ,FT_Epetra_Import_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * PermuteToLIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Import_PermuteToLIDs ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_PermuteToLIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Import_PermuteToLIDs')
    import :: c_ptr ,FT_Epetra_Import_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumRemoteIDs() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Import_NumRemoteIDs ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_NumRemoteIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Import_NumRemoteIDs')
    import :: c_int ,FT_Epetra_Import_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * RemoteLIDs() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Import_RemoteLIDs ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_RemoteLIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Import_RemoteLIDs')
    import :: c_ptr ,FT_Epetra_Import_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumExportIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Import_NumExportIDs ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_NumExportIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Import_NumExportIDs')
    import :: c_int ,FT_Epetra_Import_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * ExportLIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Import_ExportLIDs ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_ExportLIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Import_ExportLIDs')
    import :: c_ptr ,FT_Epetra_Import_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int * ExportPIDs () const;
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_Import_ExportPIDs ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_ExportPIDs ( selfID ) result(that) &
        bind(C,name='Epetra_Import_ExportPIDs')
    import :: c_ptr ,FT_Epetra_Import_ID_t
    
    type(c_ptr)                                                   :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumSend() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Import_NumSend ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_NumSend ( selfID ) result(that) &
        bind(C,name='Epetra_Import_NumSend')
    import :: c_int ,FT_Epetra_Import_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int NumRecv() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_Import_NumRecv ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_NumRecv ( selfID ) result(that) &
        bind(C,name='Epetra_Import_NumRecv')
    import :: c_int ,FT_Epetra_Import_ID_t
    
    integer(c_int)                                                :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap & SourceMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_Import_SourceMap ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_SourceMap ( selfID ) result(that) &
        bind(C,name='Epetra_Import_SourceMap')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_Import_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const Epetra_BlockMap & TargetMap() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_BlockMap_ID_t Epetra_Import_TargetMap ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_TargetMap ( selfID ) result(that) &
        bind(C,name='Epetra_Import_TargetMap')
    import :: FT_Epetra_BlockMap_ID_t ,FT_Epetra_Import_ID_t
    
    type(FT_Epetra_BlockMap_ID_t)                                  :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Distributor & Distributor() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Distributor_ID_t Epetra_Import_Distributor ( CT_Epetra_Import_ID_t selfID );

  function Epetra_Import_Distributor ( selfID ) result(that) &
        bind(C,name='Epetra_Import_Distributor')
    import :: FT_Epetra_Distributor_ID_t ,FT_Epetra_Import_ID_t
    
    type(FT_Epetra_Distributor_ID_t)                                  :: that
    type(FT_Epetra_Import_ID_t) ,intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_Time interface
!! @{

  ! _________________ Epetra_Time interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_Time_ID_t Epetra_Time_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_Time_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_Time_Degeneralize')
    import :: FT_Epetra_Time_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_Time_ID_t)                                     :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_Time_Generalize ( CT_Epetra_Time_ID_t id );

  function Epetra_Time_Generalize ( id ) result(that) bind(C,name='Epetra_Time_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_Time_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_Time_ID_t)   ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Time(const Epetra_Comm & Comm);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Time_ID_t Epetra_Time_Create ( CT_Epetra_Comm_ID_t CommID );

  function Epetra_Time_Create ( CommID ) result(that) bind(C,name='Epetra_Time_Create')
    import :: FT_Epetra_Time_ID_t ,FT_Epetra_Comm_ID_t
    
    type(FT_Epetra_Time_ID_t)                                     :: that
    type(FT_Epetra_Comm_ID_t)   ,intent(in)   ,value              :: CommID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Time(const Epetra_Time& Time);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Time_ID_t Epetra_Time_Duplicate ( CT_Epetra_Time_ID_t TimeID );

  function Epetra_Time_Duplicate ( TimeID ) result(that) &
        bind(C,name='Epetra_Time_Duplicate')
    import :: FT_Epetra_Time_ID_t
    
    type(FT_Epetra_Time_ID_t)                                     :: that
    type(FT_Epetra_Time_ID_t)   ,intent(in)   ,value              :: TimeID
  end function


  !> <BR> Original C++ prototype:
  !! double WallTime(void) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_Time_WallTime ( CT_Epetra_Time_ID_t selfID );

  function Epetra_Time_WallTime ( selfID ) result(that) bind(C,name='Epetra_Time_WallTime')
    import :: c_double ,FT_Epetra_Time_ID_t
    
    real(c_double)                                                :: that
    type(FT_Epetra_Time_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void ResetStartTime(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Time_ResetStartTime ( CT_Epetra_Time_ID_t selfID );

  subroutine Epetra_Time_ResetStartTime ( selfID ) &
        bind(C,name='Epetra_Time_ResetStartTime')
    import :: FT_Epetra_Time_ID_t
    
    type(FT_Epetra_Time_ID_t)   ,intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! double ElapsedTime(void) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_Time_ElapsedTime ( CT_Epetra_Time_ID_t selfID );

  function Epetra_Time_ElapsedTime ( selfID ) result(that) &
        bind(C,name='Epetra_Time_ElapsedTime')
    import :: c_double ,FT_Epetra_Time_ID_t
    
    real(c_double)                                                :: that
    type(FT_Epetra_Time_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_Time(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Time_Destroy ( CT_Epetra_Time_ID_t * selfID );

  subroutine Epetra_Time_Destroy ( selfID ) bind(C,name='Epetra_Time_Destroy')
    import :: FT_Epetra_Time_ID_t
    
    type(FT_Epetra_Time_ID_t)                                     :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! Epetra_Time& operator=(const Epetra_Time& src);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_Time_Assign ( CT_Epetra_Time_ID_t selfID, CT_Epetra_Time_ID_t srcID );

  subroutine Epetra_Time_Assign ( selfID, srcID ) bind(C,name='Epetra_Time_Assign')
    import :: FT_Epetra_Time_ID_t
    
    type(FT_Epetra_Time_ID_t)   ,intent(in)   ,value              :: selfID
    type(FT_Epetra_Time_ID_t)   ,intent(in)   ,value              :: srcID
  end subroutine


!> @}


!> @name Epetra_JadMatrix interface
!! @{

  ! _________________ Epetra_JadMatrix interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_JadMatrix_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_JadMatrix_Degeneralize')
    import :: FT_Epetra_JadMatrix_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_JadMatrix_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_JadMatrix_Generalize ( CT_Epetra_JadMatrix_ID_t id );

  function Epetra_JadMatrix_Generalize ( id ) result(that) &
        bind(C,name='Epetra_JadMatrix_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_JadMatrix_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_JadMatrix_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_JadMatrix(const Epetra_RowMatrix & Matrix);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_JadMatrix_ID_t Epetra_JadMatrix_Create ( CT_Epetra_RowMatrix_ID_t MatrixID );

  function Epetra_JadMatrix_Create ( MatrixID ) result(that) &
        bind(C,name='Epetra_JadMatrix_Create')
    import :: FT_Epetra_JadMatrix_ID_t ,FT_Epetra_RowMatrix_ID_t
    
    type(FT_Epetra_JadMatrix_ID_t)                                  :: that
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: MatrixID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_JadMatrix();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_JadMatrix_Destroy ( CT_Epetra_JadMatrix_ID_t * selfID );

  subroutine Epetra_JadMatrix_Destroy ( selfID ) bind(C,name='Epetra_JadMatrix_Destroy')
    import :: FT_Epetra_JadMatrix_ID_t
    
    type(FT_Epetra_JadMatrix_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int UpdateValues(const Epetra_RowMatrix & Matrix, bool CheckStructure = false);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_JadMatrix_UpdateValues ( CT_Epetra_JadMatrix_ID_t selfID, 
  !!     CT_Epetra_RowMatrix_ID_t MatrixID, boolean CheckStructure );

  function Epetra_JadMatrix_UpdateValues ( selfID, MatrixID, CheckStructure ) result(that) &
        bind(C,name='Epetra_JadMatrix_UpdateValues')
    import :: c_int ,FT_Epetra_JadMatrix_ID_t ,FT_Epetra_RowMatrix_ID_t ,FT_boolean_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_JadMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: MatrixID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: CheckStructure
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, 
  !!     int * Indices) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_JadMatrix_ExtractMyRowCopy ( CT_Epetra_JadMatrix_ID_t selfID, int MyRow, 
  !!     int Length, int * NumEntries, double * Values, int * Indices );

  function Epetra_JadMatrix_ExtractMyRowCopy ( selfID, MyRow, Length, NumEntries, Values, &
        Indices ) result(that) bind(C,name='Epetra_JadMatrix_ExtractMyRowCopy')
    import :: c_int ,FT_Epetra_JadMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_JadMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(in)   ,value              :: Length
    integer(c_int)                ,intent(inout)                    :: NumEntries
    real(c_double)                                    ,dimension(*) :: Values
    integer(c_int)                                    ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractMyEntryView(int CurEntry, double * &Value, int & RowIndex, int & ColIndex);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_JadMatrix_ExtractMyEntryView ( CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  !!     double * * Value, int * RowIndex, int * ColIndex );

  function Epetra_JadMatrix_ExtractMyEntryView ( selfID, CurEntry, Value, RowIndex, &
        ColIndex ) result(that) bind(C,name='Epetra_JadMatrix_ExtractMyEntryView')
    import :: c_int ,FT_Epetra_JadMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_JadMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: CurEntry
    real(c_double)                ,intent(inout)      ,dimension(*) :: Value
    integer(c_int)                ,intent(inout)                    :: RowIndex
    integer(c_int)                ,intent(inout)                    :: ColIndex
  end function


  !> <BR> Original C++ prototype:
  !! int ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, 
  !!     int & ColIndex) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_JadMatrix_ExtractMyEntryView_Const ( CT_Epetra_JadMatrix_ID_t selfID, int CurEntry, 
  !!     double const ** Value, int * RowIndex, int * ColIndex );

  function Epetra_JadMatrix_ExtractMyEntryView_Const ( selfID, CurEntry, Value, RowIndex, &
        ColIndex ) result(that) bind(C,name='Epetra_JadMatrix_ExtractMyEntryView_Const')
    import :: c_int ,FT_Epetra_JadMatrix_ID_t ,c_double
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_JadMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: CurEntry
    real(c_double)                ,intent(in)         ,dimension(*) :: Value
    integer(c_int)                ,intent(inout)                    :: RowIndex
    integer(c_int)                ,intent(inout)                    :: ColIndex
  end function


  !> <BR> Original C++ prototype:
  !! int NumMyRowEntries(int MyRow, int & NumEntries) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_JadMatrix_NumMyRowEntries ( CT_Epetra_JadMatrix_ID_t selfID, int MyRow, 
  !!     int * NumEntries );

  function Epetra_JadMatrix_NumMyRowEntries ( selfID, MyRow, NumEntries ) result(that) &
        bind(C,name='Epetra_JadMatrix_NumMyRowEntries')
    import :: c_int ,FT_Epetra_JadMatrix_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_JadMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                ,intent(in)   ,value              :: MyRow
    integer(c_int)                ,intent(inout)                    :: NumEntries
  end function


  !> <BR> Original C++ prototype:
  !! int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_JadMatrix_Multiply ( CT_Epetra_JadMatrix_ID_t selfID, boolean TransA, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

  function Epetra_JadMatrix_Multiply ( selfID, TransA, XID, YID ) result(that) &
        bind(C,name='Epetra_JadMatrix_Multiply')
    import :: c_int ,FT_Epetra_JadMatrix_ID_t ,FT_boolean_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_JadMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: TransA
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
  !!     Epetra_MultiVector& Y) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_JadMatrix_Solve ( CT_Epetra_JadMatrix_ID_t selfID, boolean Upper, boolean Trans, 
  !!     boolean UnitDiagonal, CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

  function Epetra_JadMatrix_Solve ( selfID, Upper, Trans, UnitDiagonal, XID, YID ) result(that) &
        bind(C,name='Epetra_JadMatrix_Solve')
    import :: c_int ,FT_Epetra_JadMatrix_ID_t ,FT_boolean_t ,FT_Epetra_MultiVector_ID_t
    
    integer(c_int)                                                  :: that
    type(FT_Epetra_JadMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)         ,intent(in)   ,value              :: Upper
    integer(FT_boolean_t)         ,intent(in)   ,value              :: Trans
    integer(FT_boolean_t)         ,intent(in)   ,value              :: UnitDiagonal
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t),intent(in)   ,value              :: YID
  end function


!> @}


!> @name Epetra_LinearProblem interface
!! @{

  ! _________________ Epetra_LinearProblem interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_LinearProblem_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_LinearProblem_Degeneralize')
    import :: FT_Epetra_LinearProblem_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t)  ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_LinearProblem_Generalize ( CT_Epetra_LinearProblem_ID_t id );

  function Epetra_LinearProblem_Generalize ( id ) result(that) &
        bind(C,name='Epetra_LinearProblem_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_LinearProblem_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                    :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_LinearProblem(void);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create (  );

  function Epetra_LinearProblem_Create (  ) result(that) &
        bind(C,name='Epetra_LinearProblem_Create')
    import :: FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t)                                  :: that
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_LinearProblem(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromMatrix ( CT_Epetra_RowMatrix_ID_t AID, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID );

  function Epetra_LinearProblem_Create_FromMatrix ( AID, XID, BID ) result(that) &
        bind(C,name='Epetra_LinearProblem_Create_FromMatrix')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_Epetra_RowMatrix_ID_t , &
          FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t)                                  :: that
    type(FT_Epetra_RowMatrix_ID_t)    ,intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t)  ,intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t)  ,intent(in)   ,value              :: BID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_LinearProblem(Epetra_Operator * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromOperator ( CT_Epetra_Operator_ID_t AID, 
  !!     CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t BID );

  function Epetra_LinearProblem_Create_FromOperator ( AID, XID, BID ) result(that) &
        bind(C,name='Epetra_LinearProblem_Create_FromOperator')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_Epetra_Operator_ID_t , &
          FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t)                                  :: that
    type(FT_Epetra_Operator_ID_t)     ,intent(in)   ,value              :: AID
    type(FT_Epetra_MultiVector_ID_t)  ,intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t)  ,intent(in)   ,value              :: BID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_LinearProblem(const Epetra_LinearProblem& Problem);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Duplicate ( CT_Epetra_LinearProblem_ID_t ProblemID );

  function Epetra_LinearProblem_Duplicate ( ProblemID ) result(that) &
        bind(C,name='Epetra_LinearProblem_Duplicate')
    import :: FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t)                                  :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: ProblemID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_LinearProblem(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LinearProblem_Destroy ( CT_Epetra_LinearProblem_ID_t * selfID );

  subroutine Epetra_LinearProblem_Destroy ( selfID ) &
        bind(C,name='Epetra_LinearProblem_Destroy')
    import :: FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int CheckInput() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_LinearProblem_CheckInput ( CT_Epetra_LinearProblem_ID_t selfID );

  function Epetra_LinearProblem_CheckInput ( selfID ) result(that) &
        bind(C,name='Epetra_LinearProblem_CheckInput')
    import :: c_int ,FT_Epetra_LinearProblem_ID_t
    
    integer(c_int)                                                      :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void AssertSymmetric();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LinearProblem_AssertSymmetric ( CT_Epetra_LinearProblem_ID_t selfID );

  subroutine Epetra_LinearProblem_AssertSymmetric ( selfID ) &
        bind(C,name='Epetra_LinearProblem_AssertSymmetric')
    import :: FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SetPDL(ProblemDifficultyLevel PDL);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LinearProblem_SetPDL ( CT_Epetra_LinearProblem_ID_t selfID, 
  !!     CT_ProblemDifficultyLevel_E_t PDL );

  subroutine Epetra_LinearProblem_SetPDL ( selfID, PDL ) &
        bind(C,name='Epetra_LinearProblem_SetPDL')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_ProblemDifficultyLevel_E_t
    
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
    integer(FT_ProblemDifficultyLevel_E_t),intent(in)   ,value              :: PDL
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SetOperator(Epetra_RowMatrix * A);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LinearProblem_SetOperator_Matrix ( CT_Epetra_LinearProblem_ID_t selfID, 
  !!     CT_Epetra_RowMatrix_ID_t AID );

  subroutine Epetra_LinearProblem_SetOperator_Matrix ( selfID, AID ) &
        bind(C,name='Epetra_LinearProblem_SetOperator_Matrix')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_Epetra_RowMatrix_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_RowMatrix_ID_t)    ,intent(in)   ,value              :: AID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SetOperator(Epetra_Operator * A);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LinearProblem_SetOperator ( CT_Epetra_LinearProblem_ID_t selfID, 
  !!     CT_Epetra_Operator_ID_t AID );

  subroutine Epetra_LinearProblem_SetOperator ( selfID, AID ) &
        bind(C,name='Epetra_LinearProblem_SetOperator')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_Epetra_Operator_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Operator_ID_t)     ,intent(in)   ,value              :: AID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SetLHS(Epetra_MultiVector * X);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LinearProblem_SetLHS ( CT_Epetra_LinearProblem_ID_t selfID, 
  !!     CT_Epetra_MultiVector_ID_t XID );

  subroutine Epetra_LinearProblem_SetLHS ( selfID, XID ) &
        bind(C,name='Epetra_LinearProblem_SetLHS')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t)  ,intent(in)   ,value              :: XID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SetRHS(Epetra_MultiVector * B);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LinearProblem_SetRHS ( CT_Epetra_LinearProblem_ID_t selfID, 
  !!     CT_Epetra_MultiVector_ID_t BID );

  subroutine Epetra_LinearProblem_SetRHS ( selfID, BID ) &
        bind(C,name='Epetra_LinearProblem_SetRHS')
    import :: FT_Epetra_LinearProblem_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t)  ,intent(in)   ,value              :: BID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int LeftScale(const Epetra_Vector & D);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_LinearProblem_LeftScale ( CT_Epetra_LinearProblem_ID_t selfID, 
  !!     CT_Epetra_Vector_ID_t DID );

  function Epetra_LinearProblem_LeftScale ( selfID, DID ) result(that) &
        bind(C,name='Epetra_LinearProblem_LeftScale')
    import :: c_int ,FT_Epetra_LinearProblem_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                      :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)       ,intent(in)   ,value              :: DID
  end function


  !> <BR> Original C++ prototype:
  !! int RightScale(const Epetra_Vector & D);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_LinearProblem_RightScale ( CT_Epetra_LinearProblem_ID_t selfID, 
  !!     CT_Epetra_Vector_ID_t DID );

  function Epetra_LinearProblem_RightScale ( selfID, DID ) result(that) &
        bind(C,name='Epetra_LinearProblem_RightScale')
    import :: c_int ,FT_Epetra_LinearProblem_ID_t ,FT_Epetra_Vector_ID_t
    
    integer(c_int)                                                      :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Vector_ID_t)       ,intent(in)   ,value              :: DID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_Operator * GetOperator() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_Operator_ID_t Epetra_LinearProblem_GetOperator ( CT_Epetra_LinearProblem_ID_t selfID );

  function Epetra_LinearProblem_GetOperator ( selfID ) result(that) &
        bind(C,name='Epetra_LinearProblem_GetOperator')
    import :: FT_Epetra_Operator_ID_t ,FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_Operator_ID_t)                                       :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_RowMatrix * GetMatrix() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_RowMatrix_ID_t Epetra_LinearProblem_GetMatrix ( CT_Epetra_LinearProblem_ID_t selfID );

  function Epetra_LinearProblem_GetMatrix ( selfID ) result(that) &
        bind(C,name='Epetra_LinearProblem_GetMatrix')
    import :: FT_Epetra_RowMatrix_ID_t ,FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_RowMatrix_ID_t)                                      :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector * GetLHS() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetLHS ( CT_Epetra_LinearProblem_ID_t selfID );

  function Epetra_LinearProblem_GetLHS ( selfID ) result(that) &
        bind(C,name='Epetra_LinearProblem_GetLHS')
    import :: FT_Epetra_MultiVector_ID_t ,FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_MultiVector_ID_t)                                    :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_MultiVector * GetRHS() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_MultiVector_ID_t Epetra_LinearProblem_GetRHS ( CT_Epetra_LinearProblem_ID_t selfID );

  function Epetra_LinearProblem_GetRHS ( selfID ) result(that) &
        bind(C,name='Epetra_LinearProblem_GetRHS')
    import :: FT_Epetra_MultiVector_ID_t ,FT_Epetra_LinearProblem_ID_t
    
    type(FT_Epetra_MultiVector_ID_t)                                    :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! ProblemDifficultyLevel GetPDL() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_ProblemDifficultyLevel_E_t Epetra_LinearProblem_GetPDL ( CT_Epetra_LinearProblem_ID_t selfID );

  function Epetra_LinearProblem_GetPDL ( selfID ) result(that) &
        bind(C,name='Epetra_LinearProblem_GetPDL')
    import :: FT_ProblemDifficultyLevel_E_t ,FT_Epetra_LinearProblem_ID_t
    
    integer(FT_ProblemDifficultyLevel_E_t)                                  :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool IsOperatorSymmetric() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_LinearProblem_IsOperatorSymmetric ( CT_Epetra_LinearProblem_ID_t selfID );

  function Epetra_LinearProblem_IsOperatorSymmetric ( selfID ) result(that) &
        bind(C,name='Epetra_LinearProblem_IsOperatorSymmetric')
    import :: FT_boolean_t ,FT_Epetra_LinearProblem_ID_t
    
    integer(FT_boolean_t)                                               :: that
    type(FT_Epetra_LinearProblem_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_LAPACK interface
!! @{

  ! _________________ Epetra_LAPACK interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_LAPACK_ID_t Epetra_LAPACK_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_LAPACK_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_LAPACK_Degeneralize')
    import :: FT_Epetra_LAPACK_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_LAPACK_ID_t)                                   :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_LAPACK_Generalize ( CT_Epetra_LAPACK_ID_t id );

  function Epetra_LAPACK_Generalize ( id ) result(that) &
        bind(C,name='Epetra_LAPACK_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_LAPACK_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_LAPACK(void);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_LAPACK_ID_t Epetra_LAPACK_Create (  );

  function Epetra_LAPACK_Create (  ) result(that) bind(C,name='Epetra_LAPACK_Create')
    import :: FT_Epetra_LAPACK_ID_t
    
    type(FT_Epetra_LAPACK_ID_t)                                   :: that
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_LAPACK(const Epetra_LAPACK& LAPACK);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_LAPACK_ID_t Epetra_LAPACK_Duplicate ( CT_Epetra_LAPACK_ID_t LAPACKID );

  function Epetra_LAPACK_Duplicate ( LAPACKID ) result(that) &
        bind(C,name='Epetra_LAPACK_Duplicate')
    import :: FT_Epetra_LAPACK_ID_t
    
    type(FT_Epetra_LAPACK_ID_t)                                   :: that
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: LAPACKID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_LAPACK(void);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_Destroy ( CT_Epetra_LAPACK_ID_t * selfID );

  subroutine Epetra_LAPACK_Destroy ( selfID ) bind(C,name='Epetra_LAPACK_Destroy')
    import :: FT_Epetra_LAPACK_ID_t
    
    type(FT_Epetra_LAPACK_ID_t)                                   :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POTRF( const char UPLO, const int N, float * A, const int LDA, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POTRF_float ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     float * A, const int LDA, int * INFO );

  subroutine Epetra_LAPACK_POTRF_float ( selfID, UPLO, N, A, LDA, INFO ) &
        bind(C,name='Epetra_LAPACK_POTRF_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POTRF( const char UPLO, const int N, double * A, const int LDA, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POTRF_double ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     double * A, const int LDA, int * INFO );

  subroutine Epetra_LAPACK_POTRF_double ( selfID, UPLO, N, A, LDA, INFO ) &
        bind(C,name='Epetra_LAPACK_POTRF_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POTRS( const char UPLO, const int N, const int NRHS, const float * A, const int LDA, 
  !!     float * X, const int LDX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POTRS_float ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     const int NRHS, const float * A, const int LDA, float * X, const int LDX, int * INFO );

  subroutine Epetra_LAPACK_POTRS_float ( selfID, UPLO, N, NRHS, A, LDA, X, LDX, INFO ) &
        bind(C,name='Epetra_LAPACK_POTRS_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POTRS( const char UPLO, const int N, const int NRHS, const double * A, const int LDA, 
  !!     double * X, const int LDX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POTRS_double ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     const int NRHS, const double * A, const int LDA, double * X, const int LDX, int * INFO );

  subroutine Epetra_LAPACK_POTRS_double ( selfID, UPLO, N, NRHS, A, LDA, X, LDX, INFO ) &
        bind(C,name='Epetra_LAPACK_POTRS_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POTRI( const char UPLO, const int N, float * A, const int LDA, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POTRI_float ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     float * A, const int LDA, int * INFO );

  subroutine Epetra_LAPACK_POTRI_float ( selfID, UPLO, N, A, LDA, INFO ) &
        bind(C,name='Epetra_LAPACK_POTRI_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POTRI( const char UPLO, const int N, double * A, const int LDA, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POTRI_double ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     double * A, const int LDA, int * INFO );

  subroutine Epetra_LAPACK_POTRI_double ( selfID, UPLO, N, A, LDA, INFO ) &
        bind(C,name='Epetra_LAPACK_POTRI_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POCON( const char UPLO, const int N, const float * A, const int LDA, const float ANORM, 
  !!     float * RCOND, float * WORK, int * IWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POCON_float ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     const float * A, const int LDA, const float ANORM, float * RCOND, float * WORK, 
  !!     int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_POCON_float ( selfID, UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK, &
        INFO ) bind(C,name='Epetra_LAPACK_POCON_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)   ,value              :: ANORM
    real(c_float)                                   ,dimension(*) :: RCOND
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POCON( const char UPLO, const int N, const double * A, const int LDA, const double ANORM, 
  !!     double * RCOND, double * WORK, int * IWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POCON_double ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     const double * A, const int LDA, const double ANORM, double * RCOND, double * WORK, 
  !!     int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_POCON_double ( selfID, UPLO, N, A, LDA, ANORM, RCOND, WORK, &
        IWORK, INFO ) bind(C,name='Epetra_LAPACK_POCON_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)   ,value              :: ANORM
    real(c_double)                                  ,dimension(*) :: RCOND
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POSV( const char UPLO, const int N, const int NRHS, float * A, const int LDA, float * X, 
  !!     const int LDX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POSV_float ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     const int NRHS, float * A, const int LDA, float * X, const int LDX, int * INFO );

  subroutine Epetra_LAPACK_POSV_float ( selfID, UPLO, N, NRHS, A, LDA, X, LDX, INFO ) &
        bind(C,name='Epetra_LAPACK_POSV_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POSV( const char UPLO, const int N, const int NRHS, double * A, const int LDA, 
  !!     double * X, const int LDX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POSV_double ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     const int NRHS, double * A, const int LDA, double * X, const int LDX, int * INFO );

  subroutine Epetra_LAPACK_POSV_double ( selfID, UPLO, N, NRHS, A, LDA, X, LDX, INFO ) &
        bind(C,name='Epetra_LAPACK_POSV_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POEQU(const int N, const float * A, const int LDA, float * S, float * SCOND, 
  !!     float * AMAX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POEQU_float ( CT_Epetra_LAPACK_ID_t selfID, const int N, const float * A, 
  !!     const int LDA, float * S, float * SCOND, float * AMAX, int * INFO );

  subroutine Epetra_LAPACK_POEQU_float ( selfID, N, A, LDA, S, SCOND, AMAX, INFO ) &
        bind(C,name='Epetra_LAPACK_POEQU_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: S
    real(c_float)                                   ,dimension(*) :: SCOND
    real(c_float)                                   ,dimension(*) :: AMAX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POEQU(const int N, const double * A, const int LDA, double * S, double * SCOND, 
  !!     double * AMAX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POEQU_double ( CT_Epetra_LAPACK_ID_t selfID, const int N, const double * A, 
  !!     const int LDA, double * S, double * SCOND, double * AMAX, int * INFO );

  subroutine Epetra_LAPACK_POEQU_double ( selfID, N, A, LDA, S, SCOND, AMAX, INFO ) &
        bind(C,name='Epetra_LAPACK_POEQU_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: S
    real(c_double)                                  ,dimension(*) :: SCOND
    real(c_double)                                  ,dimension(*) :: AMAX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void PORFS(const char UPLO, const int N, const int NRHS, const float * A, const int LDA, 
  !!     const float * AF, const int LDAF, const float * B, const int LDB, float * X, 
  !!     const int LDX, float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_PORFS_float ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     const int NRHS, const float * A, const int LDA, const float * AF, const int LDAF, 
  !!     const float * B, const int LDB, float * X, const int LDX, float * FERR, float * BERR, 
  !!     float * WORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_PORFS_float ( selfID, UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, &
        LDX, FERR, BERR, WORK, IWORK, INFO ) bind(C,name='Epetra_LAPACK_PORFS_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)         ,dimension(*) :: AF
    integer(c_int)              ,intent(in)   ,value              :: LDAF
    real(c_float)               ,intent(in)         ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)                                   ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    real(c_float)                                   ,dimension(*) :: FERR
    real(c_float)                                   ,dimension(*) :: BERR
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void PORFS(const char UPLO, const int N, const int NRHS, const double * A, const int LDA, 
  !!     const double * AF, const int LDAF, const double * B, const int LDB, double * X, 
  !!     const int LDX, double * FERR, double * BERR, double * WORK, int * IWORK, 
  !!     int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_PORFS_double ( CT_Epetra_LAPACK_ID_t selfID, const char UPLO, const int N, 
  !!     const int NRHS, const double * A, const int LDA, const double * AF, const int LDAF, 
  !!     const double * B, const int LDB, double * X, const int LDX, double * FERR, double * BERR, 
  !!     double * WORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_PORFS_double ( selfID, UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, &
        X, LDX, FERR, BERR, WORK, IWORK, INFO ) bind(C,name='Epetra_LAPACK_PORFS_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)         ,dimension(*) :: AF
    integer(c_int)              ,intent(in)   ,value              :: LDAF
    real(c_double)              ,intent(in)         ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)                                  ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    real(c_double)                                  ,dimension(*) :: FERR
    real(c_double)                                  ,dimension(*) :: BERR
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POSVX(const char FACT, const char UPLO, const int N, const int NRHS, float * A, 
  !!     const int LDA, float * AF, const int LDAF, const char EQUED, float * S, float * B, 
  !!     const int LDB, float * X, const int LDX, float * RCOND, float * FERR, float * BERR, 
  !!     float * WORK, int * IWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POSVX_float ( CT_Epetra_LAPACK_ID_t selfID, const char FACT, 
  !!     const char UPLO, const int N, const int NRHS, float * A, const int LDA, float * AF, 
  !!     const int LDAF, const char EQUED, float * S, float * B, const int LDB, float * X, 
  !!     const int LDX, float * RCOND, float * FERR, float * BERR, float * WORK, int * IWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_POSVX_float ( selfID, FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, &
        EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_POSVX_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: FACT
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: AF
    integer(c_int)              ,intent(in)   ,value              :: LDAF
    character(kind=c_char)      ,intent(in)   ,value              :: EQUED
    real(c_float)                                   ,dimension(*) :: S
    real(c_float)                                   ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)                                   ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    real(c_float)                                   ,dimension(*) :: RCOND
    real(c_float)                                   ,dimension(*) :: FERR
    real(c_float)                                   ,dimension(*) :: BERR
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void POSVX(const char FACT, const char UPLO, const int N, const int NRHS, double * A, 
  !!     const int LDA, double * AF, const int LDAF, const char EQUED, double * S, double * B, 
  !!     const int LDB, double * X, const int LDX, double * RCOND, double * FERR, double * BERR, 
  !!     double * WORK, int * IWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_POSVX_double ( CT_Epetra_LAPACK_ID_t selfID, const char FACT, 
  !!     const char UPLO, const int N, const int NRHS, double * A, const int LDA, double * AF, 
  !!     const int LDAF, const char EQUED, double * S, double * B, const int LDB, double * X, 
  !!     const int LDX, double * RCOND, double * FERR, double * BERR, double * WORK, int * IWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_POSVX_double ( selfID, FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, &
        EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_POSVX_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: FACT
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: AF
    integer(c_int)              ,intent(in)   ,value              :: LDAF
    character(kind=c_char)      ,intent(in)   ,value              :: EQUED
    real(c_double)                                  ,dimension(*) :: S
    real(c_double)                                  ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)                                  ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    real(c_double)                                  ,dimension(*) :: RCOND
    real(c_double)                                  ,dimension(*) :: FERR
    real(c_double)                                  ,dimension(*) :: BERR
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GELS( const char TRANS, const int M, const int N, const int NRHS, double* A, 
  !!     const int LDA, double* B, const int LDB, double* WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GELS_double ( CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int M, 
  !!     const int N, const int NRHS, double * A, const int LDA, double * B, const int LDB, 
  !!     double * WORK, const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_GELS_double ( selfID, TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, &
        LWORK, INFO ) bind(C,name='Epetra_LAPACK_GELS_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GETRF( const int M, const int N, float * A, const int LDA, int * IPIV, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GETRF_float ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     float * A, const int LDA, int * IPIV, int * INFO );

  subroutine Epetra_LAPACK_GETRF_float ( selfID, M, N, A, LDA, IPIV, INFO ) &
        bind(C,name='Epetra_LAPACK_GETRF_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: IPIV
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GETRF( const int M, const int N, double * A, const int LDA, int * IPIV, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GETRF_double ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     double * A, const int LDA, int * IPIV, int * INFO );

  subroutine Epetra_LAPACK_GETRF_double ( selfID, M, N, A, LDA, IPIV, INFO ) &
        bind(C,name='Epetra_LAPACK_GETRF_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: IPIV
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEQRF( const int M, const int N, float * A, const int LDA, float * TAU, float * WORK, 
  !!     const int lwork, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEQRF_float ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     float * A, const int LDA, float * TAU, float * WORK, const int lwork, int * INFO );

  subroutine Epetra_LAPACK_GEQRF_float ( selfID, M, N, A, LDA, TAU, WORK, lwork, INFO ) &
        bind(C,name='Epetra_LAPACK_GEQRF_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: TAU
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: lwork
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEQRF( const int M, const int N, double * A, const int LDA, double * TAU, double * WORK, 
  !!     const int lwork, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEQRF_double ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     double * A, const int LDA, double * TAU, double * WORK, const int lwork, int * INFO );

  subroutine Epetra_LAPACK_GEQRF_double ( selfID, M, N, A, LDA, TAU, WORK, lwork, INFO ) &
        bind(C,name='Epetra_LAPACK_GEQRF_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: TAU
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: lwork
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GETRS( const char TRANS, const int N, const int NRHS, const float * A, const int LDA, 
  !!     const int * IPIV, float * X, const int LDX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GETRS_float ( CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int N, 
  !!     const int NRHS, const float * A, const int LDA, const int * IPIV, float * X, 
  !!     const int LDX, int * INFO );

  subroutine Epetra_LAPACK_GETRS_float ( selfID, TRANS, N, NRHS, A, LDA, IPIV, X, LDX, INFO ) &
        bind(C,name='Epetra_LAPACK_GETRS_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)              ,intent(in)         ,dimension(*) :: IPIV
    real(c_float)                                   ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GETRS( const char TRANS, const int N, const int NRHS, const double * A, const int LDA, 
  !!     const int * IPIV, double * X, const int LDX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GETRS_double ( CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int N, 
  !!     const int NRHS, const double * A, const int LDA, const int * IPIV, double * X, 
  !!     const int LDX, int * INFO );

  subroutine Epetra_LAPACK_GETRS_double ( selfID, TRANS, N, NRHS, A, LDA, IPIV, X, LDX, &
        INFO ) bind(C,name='Epetra_LAPACK_GETRS_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)              ,intent(in)         ,dimension(*) :: IPIV
    real(c_double)                                  ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GETRI( const int N, float * A, const int LDA, int * IPIV, float * WORK, 
  !!     const int * LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GETRI_float ( CT_Epetra_LAPACK_ID_t selfID, const int N, float * A, 
  !!     const int LDA, int * IPIV, float * WORK, const int * LWORK, int * INFO );

  subroutine Epetra_LAPACK_GETRI_float ( selfID, N, A, LDA, IPIV, WORK, LWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GETRI_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: IPIV
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)         ,dimension(*) :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GETRI( const int N, double * A, const int LDA, int * IPIV, double * WORK, 
  !!     const int * LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GETRI_double ( CT_Epetra_LAPACK_ID_t selfID, const int N, double * A, 
  !!     const int LDA, int * IPIV, double * WORK, const int * LWORK, int * INFO );

  subroutine Epetra_LAPACK_GETRI_double ( selfID, N, A, LDA, IPIV, WORK, LWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GETRI_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: IPIV
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)         ,dimension(*) :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GECON( const char NORM, const int N, const float * A, const int LDA, const float ANORM, 
  !!     float * RCOND, float * WORK, int * IWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GECON_float ( CT_Epetra_LAPACK_ID_t selfID, const char NORM, const int N, 
  !!     const float * A, const int LDA, const float ANORM, float * RCOND, float * WORK, 
  !!     int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GECON_float ( selfID, NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, &
        INFO ) bind(C,name='Epetra_LAPACK_GECON_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: NORM
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)   ,value              :: ANORM
    real(c_float)                                   ,dimension(*) :: RCOND
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GECON( const char NORM, const int N, const double * A, const int LDA, const double ANORM, 
  !!     double * RCOND, double * WORK, int * IWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GECON_double ( CT_Epetra_LAPACK_ID_t selfID, const char NORM, const int N, 
  !!     const double * A, const int LDA, const double ANORM, double * RCOND, double * WORK, 
  !!     int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GECON_double ( selfID, NORM, N, A, LDA, ANORM, RCOND, WORK, &
        IWORK, INFO ) bind(C,name='Epetra_LAPACK_GECON_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: NORM
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)   ,value              :: ANORM
    real(c_double)                                  ,dimension(*) :: RCOND
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GESV( const int N, const int NRHS, float * A, const int LDA, int * IPIV, float * X, 
  !!     const int LDX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GESV_float ( CT_Epetra_LAPACK_ID_t selfID, const int N, const int NRHS, 
  !!     float * A, const int LDA, int * IPIV, float * X, const int LDX, int * INFO );

  subroutine Epetra_LAPACK_GESV_float ( selfID, N, NRHS, A, LDA, IPIV, X, LDX, INFO ) &
        bind(C,name='Epetra_LAPACK_GESV_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: IPIV
    real(c_float)                                   ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GESV( const int N, const int NRHS, double * A, const int LDA, int * IPIV, double * X, 
  !!     const int LDX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GESV_double ( CT_Epetra_LAPACK_ID_t selfID, const int N, const int NRHS, 
  !!     double * A, const int LDA, int * IPIV, double * X, const int LDX, int * INFO );

  subroutine Epetra_LAPACK_GESV_double ( selfID, N, NRHS, A, LDA, IPIV, X, LDX, INFO ) &
        bind(C,name='Epetra_LAPACK_GESV_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    integer(c_int)                                  ,dimension(*) :: IPIV
    real(c_double)                                  ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEEQU(const int M, const int N, const float * A, const int LDA, float * R, float * C, 
  !!     float * ROWCND, float * COLCND, float * AMAX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEEQU_float ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     const float * A, const int LDA, float * R, float * C, float * ROWCND, float * COLCND, 
  !!     float * AMAX, int * INFO );

  subroutine Epetra_LAPACK_GEEQU_float ( selfID, M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
        INFO ) bind(C,name='Epetra_LAPACK_GEEQU_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: R
    real(c_float)                                   ,dimension(*) :: C
    real(c_float)                                   ,dimension(*) :: ROWCND
    real(c_float)                                   ,dimension(*) :: COLCND
    real(c_float)                                   ,dimension(*) :: AMAX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEEQU(const int M, const int N, const double * A, const int LDA, double * R, double * C, 
  !!     double * ROWCND, double * COLCND, double * AMAX, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEEQU_double ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     const double * A, const int LDA, double * R, double * C, double * ROWCND, 
  !!     double * COLCND, double * AMAX, int * INFO );

  subroutine Epetra_LAPACK_GEEQU_double ( selfID, M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, &
        INFO ) bind(C,name='Epetra_LAPACK_GEEQU_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: R
    real(c_double)                                  ,dimension(*) :: C
    real(c_double)                                  ,dimension(*) :: ROWCND
    real(c_double)                                  ,dimension(*) :: COLCND
    real(c_double)                                  ,dimension(*) :: AMAX
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GERFS(const char TRANS, const int N, const int NRHS, const float * A, const int LDA, 
  !!     const float * AF, const int LDAF, const int * IPIV, const float * B, const int LDB, 
  !!     float * X, const int LDX, float * FERR, float * BERR, float * WORK, int * IWORK, 
  !!     int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GERFS_float ( CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int N, 
  !!     const int NRHS, const float * A, const int LDA, const float * AF, const int LDAF, 
  !!     const int * IPIV, const float * B, const int LDB, float * X, const int LDX, float * FERR, 
  !!     float * BERR, float * WORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GERFS_float ( selfID, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, &
        LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GERFS_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)         ,dimension(*) :: AF
    integer(c_int)              ,intent(in)   ,value              :: LDAF
    integer(c_int)              ,intent(in)         ,dimension(*) :: IPIV
    real(c_float)               ,intent(in)         ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)                                   ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    real(c_float)                                   ,dimension(*) :: FERR
    real(c_float)                                   ,dimension(*) :: BERR
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GERFS(const char TRANS, const int N, const int NRHS, const double * A, const int LDA, 
  !!     const double * AF, const int LDAF, const int * IPIV, const double * B, const int LDB, 
  !!     double * X, const int LDX, double * FERR, double * BERR, double * WORK, int * IWORK, 
  !!     int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GERFS_double ( CT_Epetra_LAPACK_ID_t selfID, const char TRANS, const int N, 
  !!     const int NRHS, const double * A, const int LDA, const double * AF, const int LDAF, 
  !!     const int * IPIV, const double * B, const int LDB, double * X, const int LDX, 
  !!     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GERFS_double ( selfID, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, &
        LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GERFS_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)         ,dimension(*) :: AF
    integer(c_int)              ,intent(in)   ,value              :: LDAF
    integer(c_int)              ,intent(in)         ,dimension(*) :: IPIV
    real(c_double)              ,intent(in)         ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)                                  ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    real(c_double)                                  ,dimension(*) :: FERR
    real(c_double)                                  ,dimension(*) :: BERR
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GESVX(const char FACT, const char TRANS, const int N, const int NRHS, float * A, 
  !!     const int LDA, float * AF, const int LDAF, int * IPIV, const char EQUED, float * R, 
  !!     float * C, float * B, const int LDB, float * X, const int LDX, float * RCOND, 
  !!     float * FERR, float * BERR, float * WORK, int * IWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GESVX_float ( CT_Epetra_LAPACK_ID_t selfID, const char FACT, 
  !!     const char TRANS, const int N, const int NRHS, float * A, const int LDA, float * AF, 
  !!     const int LDAF, int * IPIV, const char EQUED, float * R, float * C, float * B, 
  !!     const int LDB, float * X, const int LDX, float * RCOND, float * FERR, float * BERR, 
  !!     float * WORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GESVX_float ( selfID, FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, &
        IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GESVX_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: FACT
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: AF
    integer(c_int)              ,intent(in)   ,value              :: LDAF
    integer(c_int)                                  ,dimension(*) :: IPIV
    character(kind=c_char)      ,intent(in)   ,value              :: EQUED
    real(c_float)                                   ,dimension(*) :: R
    real(c_float)                                   ,dimension(*) :: C
    real(c_float)                                   ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)                                   ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    real(c_float)                                   ,dimension(*) :: RCOND
    real(c_float)                                   ,dimension(*) :: FERR
    real(c_float)                                   ,dimension(*) :: BERR
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GESVX(const char FACT, const char TRANS, const int N, const int NRHS, double * A, 
  !!     const int LDA, double * AF, const int LDAF, int * IPIV, const char EQUED, double * R, 
  !!     double * C, double * B, const int LDB, double * X, const int LDX, double * RCOND, 
  !!     double * FERR, double * BERR, double * WORK, int * IWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GESVX_double ( CT_Epetra_LAPACK_ID_t selfID, const char FACT, 
  !!     const char TRANS, const int N, const int NRHS, double * A, const int LDA, double * AF, 
  !!     const int LDAF, int * IPIV, const char EQUED, double * R, double * C, double * B, 
  !!     const int LDB, double * X, const int LDX, double * RCOND, double * FERR, double * BERR, 
  !!     double * WORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GESVX_double ( selfID, FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, &
        IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GESVX_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: FACT
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: NRHS
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: AF
    integer(c_int)              ,intent(in)   ,value              :: LDAF
    integer(c_int)                                  ,dimension(*) :: IPIV
    character(kind=c_char)      ,intent(in)   ,value              :: EQUED
    real(c_double)                                  ,dimension(*) :: R
    real(c_double)                                  ,dimension(*) :: C
    real(c_double)                                  ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)                                  ,dimension(*) :: X
    integer(c_int)              ,intent(in)   ,value              :: LDX
    real(c_double)                                  ,dimension(*) :: RCOND
    real(c_double)                                  ,dimension(*) :: FERR
    real(c_double)                                  ,dimension(*) :: BERR
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEHRD(const int N, const int ILO, const int IHI, float * A, const int LDA, float * TAU, 
  !!     float * WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEHRD_float ( CT_Epetra_LAPACK_ID_t selfID, const int N, const int ILO, 
  !!     const int IHI, float * A, const int LDA, float * TAU, float * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_GEHRD_float ( selfID, N, ILO, IHI, A, LDA, TAU, WORK, LWORK, &
        INFO ) bind(C,name='Epetra_LAPACK_GEHRD_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: ILO
    integer(c_int)              ,intent(in)   ,value              :: IHI
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: TAU
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEHRD(const int N, const int ILO, const int IHI, double * A, const int LDA, double * TAU, 
  !!     double * WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEHRD_double ( CT_Epetra_LAPACK_ID_t selfID, const int N, const int ILO, 
  !!     const int IHI, double * A, const int LDA, double * TAU, double * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_GEHRD_double ( selfID, N, ILO, IHI, A, LDA, TAU, WORK, LWORK, &
        INFO ) bind(C,name='Epetra_LAPACK_GEHRD_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: ILO
    integer(c_int)              ,intent(in)   ,value              :: IHI
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: TAU
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void HSEQR( const char JOB, const char COMPZ, const int N, const int ILO, const int IHI, 
  !!     float * H, const int LDH, float * WR, float * WI, float * Z, const int LDZ, float * WORK, 
  !!     const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_HSEQR_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOB, 
  !!     const char COMPZ, const int N, const int ILO, const int IHI, float * H, const int LDH, 
  !!     float * WR, float * WI, float * Z, const int LDZ, float * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_HSEQR_float ( selfID, JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, &
        LDZ, WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_HSEQR_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOB
    character(kind=c_char)      ,intent(in)   ,value              :: COMPZ
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: ILO
    integer(c_int)              ,intent(in)   ,value              :: IHI
    real(c_float)                                   ,dimension(*) :: H
    integer(c_int)              ,intent(in)   ,value              :: LDH
    real(c_float)                                   ,dimension(*) :: WR
    real(c_float)                                   ,dimension(*) :: WI
    real(c_float)                                   ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void HSEQR( const char JOB, const char COMPZ, const int N, const int ILO, const int IHI, 
  !!     double * H, const int LDH, double * WR, double * WI, double * Z, const int LDZ, 
  !!     double * WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_HSEQR_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOB, 
  !!     const char COMPZ, const int N, const int ILO, const int IHI, double * H, const int LDH, 
  !!     double * WR, double * WI, double * Z, const int LDZ, double * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_HSEQR_double ( selfID, JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, &
        Z, LDZ, WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_HSEQR_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOB
    character(kind=c_char)      ,intent(in)   ,value              :: COMPZ
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: ILO
    integer(c_int)              ,intent(in)   ,value              :: IHI
    real(c_double)                                  ,dimension(*) :: H
    integer(c_int)              ,intent(in)   ,value              :: LDH
    real(c_double)                                  ,dimension(*) :: WR
    real(c_double)                                  ,dimension(*) :: WI
    real(c_double)                                  ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void ORGQR( const int M, const int N, const int K, float * A, const int LDA, float * TAU, 
  !!     float * WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_ORGQR_float ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     const int K, float * A, const int LDA, float * TAU, float * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_ORGQR_float ( selfID, M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_ORGQR_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: K
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: TAU
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void ORGQR( const int M, const int N, const int K, double * A, const int LDA, double * TAU, 
  !!     double * WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_ORGQR_double ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     const int K, double * A, const int LDA, double * TAU, double * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_ORGQR_double ( selfID, M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_ORGQR_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: K
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: TAU
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void ORGHR( const int N, const int ILO, const int IHI, float * A, const int LDA, float * TAU, 
  !!     float * WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_ORGHR_float ( CT_Epetra_LAPACK_ID_t selfID, const int N, const int ILO, 
  !!     const int IHI, float * A, const int LDA, float * TAU, float * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_ORGHR_float ( selfID, N, ILO, IHI, A, LDA, TAU, WORK, LWORK, &
        INFO ) bind(C,name='Epetra_LAPACK_ORGHR_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: ILO
    integer(c_int)              ,intent(in)   ,value              :: IHI
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: TAU
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void ORGHR( const int N, const int ILO, const int IHI, double * A, const int LDA, 
  !!     double * TAU, double * WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_ORGHR_double ( CT_Epetra_LAPACK_ID_t selfID, const int N, const int ILO, 
  !!     const int IHI, double * A, const int LDA, double * TAU, double * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_ORGHR_double ( selfID, N, ILO, IHI, A, LDA, TAU, WORK, LWORK, &
        INFO ) bind(C,name='Epetra_LAPACK_ORGHR_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: ILO
    integer(c_int)              ,intent(in)   ,value              :: IHI
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: TAU
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void ORMHR( const char SIDE, const char TRANS, const int M, const int N, const int ILO, 
  !!     const int IHI, const float * A, const int LDA, const float * TAU, float * C, 
  !!     const int LDC, float * WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_ORMHR_float ( CT_Epetra_LAPACK_ID_t selfID, const char SIDE, 
  !!     const char TRANS, const int M, const int N, const int ILO, const int IHI, 
  !!     const float * A, const int LDA, const float * TAU, float * C, const int LDC, 
  !!     float * WORK, const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_ORMHR_float ( selfID, SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, &
        C, LDC, WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_ORMHR_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: SIDE
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: ILO
    integer(c_int)              ,intent(in)   ,value              :: IHI
    real(c_float)               ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)         ,dimension(*) :: TAU
    real(c_float)                                   ,dimension(*) :: C
    integer(c_int)              ,intent(in)   ,value              :: LDC
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void ORMHR( const char SIDE, const char TRANS, const int M, const int N, const int ILO, 
  !!     const int IHI, const double * A, const int LDA, const double * TAU, double * C, 
  !!     const int LDC, double * WORK, const int LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_ORMHR_double ( CT_Epetra_LAPACK_ID_t selfID, const char SIDE, 
  !!     const char TRANS, const int M, const int N, const int ILO, const int IHI, 
  !!     const double * A, const int LDA, const double * TAU, double * C, const int LDC, 
  !!     double * WORK, const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_ORMHR_double ( selfID, SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, &
        C, LDC, WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_ORMHR_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: SIDE
    character(kind=c_char)      ,intent(in)   ,value              :: TRANS
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: ILO
    integer(c_int)              ,intent(in)   ,value              :: IHI
    real(c_double)              ,intent(in)         ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)         ,dimension(*) :: TAU
    real(c_double)                                  ,dimension(*) :: C
    integer(c_int)              ,intent(in)   ,value              :: LDC
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void LARFT( const char DIRECT, const char STOREV, const int N, const int K, double * V, 
  !!     const int LDV, double * TAU, double * T, const int LDT) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_LARFT_float ( CT_Epetra_LAPACK_ID_t selfID, const char DIRECT, 
  !!     const char STOREV, const int N, const int K, double * V, const int LDV, double * TAU, 
  !!     double * T, const int LDT );

  subroutine Epetra_LAPACK_LARFT_float ( selfID, DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) &
        bind(C,name='Epetra_LAPACK_LARFT_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: DIRECT
    character(kind=c_char)      ,intent(in)   ,value              :: STOREV
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: K
    real(c_double)                                  ,dimension(*) :: V
    integer(c_int)              ,intent(in)   ,value              :: LDV
    real(c_double)                                  ,dimension(*) :: TAU
    real(c_double)                                  ,dimension(*) :: T
    integer(c_int)              ,intent(in)   ,value              :: LDT
  end subroutine


  !> <BR> Original C++ prototype:
  !! void LARFT( const char DIRECT, const char STOREV, const int N, const int K, float * V, 
  !!     const int LDV, float * TAU, float * T, const int LDT) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_LARFT_double ( CT_Epetra_LAPACK_ID_t selfID, const char DIRECT, 
  !!     const char STOREV, const int N, const int K, float * V, const int LDV, float * TAU, 
  !!     float * T, const int LDT );

  subroutine Epetra_LAPACK_LARFT_double ( selfID, DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) &
        bind(C,name='Epetra_LAPACK_LARFT_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: DIRECT
    character(kind=c_char)      ,intent(in)   ,value              :: STOREV
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: K
    real(c_float)                                   ,dimension(*) :: V
    integer(c_int)              ,intent(in)   ,value              :: LDV
    real(c_float)                                   ,dimension(*) :: TAU
    real(c_float)                                   ,dimension(*) :: T
    integer(c_int)              ,intent(in)   ,value              :: LDT
  end subroutine


  !> <BR> Original C++ prototype:
  !! void TREVC( const char SIDE, const char HOWMNY, int * SELECT, const int N, const float * T, 
  !!     const int LDT, float *VL, const int LDVL, float * VR, const int LDVR, const int MM, 
  !!     int * M, float * WORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_TREVC_float ( CT_Epetra_LAPACK_ID_t selfID, const char SIDE, 
  !!     const char HOWMNY, int * SELECT, const int N, const float * T, const int LDT, float * VL, 
  !!     const int LDVL, float * VR, const int LDVR, const int MM, int * M, float * WORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_TREVC_float ( selfID, SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, &
        VR, LDVR, MM, M, WORK, INFO ) bind(C,name='Epetra_LAPACK_TREVC_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: SIDE
    character(kind=c_char)      ,intent(in)   ,value              :: HOWMNY
    integer(c_int)                                  ,dimension(*) :: SELECT
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)               ,intent(in)         ,dimension(*) :: T
    integer(c_int)              ,intent(in)   ,value              :: LDT
    real(c_float)                                   ,dimension(*) :: VL
    integer(c_int)              ,intent(in)   ,value              :: LDVL
    real(c_float)                                   ,dimension(*) :: VR
    integer(c_int)              ,intent(in)   ,value              :: LDVR
    integer(c_int)              ,intent(in)   ,value              :: MM
    integer(c_int)                                  ,dimension(*) :: M
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void TREVC( const char SIDE, const char HOWMNY, int * SELECT, const int N, const double * T, 
  !!     const int LDT, double *VL, const int LDVL, double * VR, const int LDVR, const int MM, 
  !!     int *M, double * WORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_TREVC_double ( CT_Epetra_LAPACK_ID_t selfID, const char SIDE, 
  !!     const char HOWMNY, int * SELECT, const int N, const double * T, const int LDT, 
  !!     double * VL, const int LDVL, double * VR, const int LDVR, const int MM, int * M, 
  !!     double * WORK, int * INFO );

  subroutine Epetra_LAPACK_TREVC_double ( selfID, SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, &
        VR, LDVR, MM, M, WORK, INFO ) bind(C,name='Epetra_LAPACK_TREVC_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: SIDE
    character(kind=c_char)      ,intent(in)   ,value              :: HOWMNY
    integer(c_int)                                  ,dimension(*) :: SELECT
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)              ,intent(in)         ,dimension(*) :: T
    integer(c_int)              ,intent(in)   ,value              :: LDT
    real(c_double)                                  ,dimension(*) :: VL
    integer(c_int)              ,intent(in)   ,value              :: LDVL
    real(c_double)                                  ,dimension(*) :: VR
    integer(c_int)              ,intent(in)   ,value              :: LDVR
    integer(c_int)              ,intent(in)   ,value              :: MM
    integer(c_int)                                  ,dimension(*) :: M
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void TREXC( const char COMPQ, const int N, float * T, const int LDT, float * Q, const int LDQ, 
  !!     int IFST, int ILST, float * WORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_TREXC_float ( CT_Epetra_LAPACK_ID_t selfID, const char COMPQ, const int N, 
  !!     float * T, const int LDT, float * Q, const int LDQ, int IFST, int ILST, float * WORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_TREXC_float ( selfID, COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, &
        INFO ) bind(C,name='Epetra_LAPACK_TREXC_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: COMPQ
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: T
    integer(c_int)              ,intent(in)   ,value              :: LDT
    real(c_float)                                   ,dimension(*) :: Q
    integer(c_int)              ,intent(in)   ,value              :: LDQ
    integer(c_int)              ,intent(in)   ,value              :: IFST
    integer(c_int)              ,intent(in)   ,value              :: ILST
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void TREXC( const char COMPQ, const int N, double * T, const int LDT, double * Q, 
  !!     const int LDQ, int IFST, int ILST, double * WORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_TREXC_double ( CT_Epetra_LAPACK_ID_t selfID, const char COMPQ, const int N, 
  !!     double * T, const int LDT, double * Q, const int LDQ, int IFST, int ILST, double * WORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_TREXC_double ( selfID, COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, &
        WORK, INFO ) bind(C,name='Epetra_LAPACK_TREXC_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: COMPQ
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: T
    integer(c_int)              ,intent(in)   ,value              :: LDT
    real(c_double)                                  ,dimension(*) :: Q
    integer(c_int)              ,intent(in)   ,value              :: LDQ
    integer(c_int)              ,intent(in)   ,value              :: IFST
    integer(c_int)              ,intent(in)   ,value              :: ILST
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GESVD( const char JOBU, const char JOBVT, const int M, const int N, float * A, 
  !!     const int LDA, float * S, float * U, const int LDU, float * VT, const int LDVT, 
  !!     float * WORK, const int * LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GESVD_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBU, 
  !!     const char JOBVT, const int M, const int N, float * A, const int LDA, float * S, 
  !!     float * U, const int LDU, float * VT, const int LDVT, float * WORK, const int * LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_GESVD_float ( selfID, JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, &
        LDVT, WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_GESVD_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBU
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVT
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: S
    real(c_float)                                   ,dimension(*) :: U
    integer(c_int)              ,intent(in)   ,value              :: LDU
    real(c_float)                                   ,dimension(*) :: VT
    integer(c_int)              ,intent(in)   ,value              :: LDVT
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)         ,dimension(*) :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GESVD( const char JOBU, const char JOBVT, const int M, const int N, double * A, 
  !!     const int LDA, double * S, double * U, const int LDU, double * VT, const int LDVT, 
  !!     double * WORK, const int * LWORK, int * INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GESVD_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBU, 
  !!     const char JOBVT, const int M, const int N, double * A, const int LDA, double * S, 
  !!     double * U, const int LDU, double * VT, const int LDVT, double * WORK, const int * LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_GESVD_double ( selfID, JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, &
        LDVT, WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_GESVD_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBU
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVT
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: S
    real(c_double)                                  ,dimension(*) :: U
    integer(c_int)              ,intent(in)   ,value              :: LDU
    real(c_double)                                  ,dimension(*) :: VT
    integer(c_int)              ,intent(in)   ,value              :: LDVT
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)         ,dimension(*) :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GGSVD(const char JOBU, const char JOBV, const char JOBQ, const int M, const int N, 
  !!     const int P, int * K, int * L, double* A, const int LDA, double* B, const int LDB, 
  !!     double* ALPHA, double* BETA, double* U, const int LDU, double* V, const int LDV, 
  !!     double* Q, const int LDQ, double* WORK, int* IWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GGSVD_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBU, 
  !!     const char JOBV, const char JOBQ, const int M, const int N, const int P, int * K, 
  !!     int * L, double * A, const int LDA, double * B, const int LDB, double * ALPHA, 
  !!     double * BETA, double * U, const int LDU, double * V, const int LDV, double * Q, 
  !!     const int LDQ, double * WORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GGSVD_double ( selfID, JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, &
        B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, IWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GGSVD_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBU
    character(kind=c_char)      ,intent(in)   ,value              :: JOBV
    character(kind=c_char)      ,intent(in)   ,value              :: JOBQ
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: P
    integer(c_int)                                  ,dimension(*) :: K
    integer(c_int)                                  ,dimension(*) :: L
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)                                  ,dimension(*) :: ALPHA
    real(c_double)                                  ,dimension(*) :: BETA
    real(c_double)                                  ,dimension(*) :: U
    integer(c_int)              ,intent(in)   ,value              :: LDU
    real(c_double)                                  ,dimension(*) :: V
    integer(c_int)              ,intent(in)   ,value              :: LDV
    real(c_double)                                  ,dimension(*) :: Q
    integer(c_int)              ,intent(in)   ,value              :: LDQ
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GGSVD(const char JOBU, const char JOBV, const char JOBQ, const int M, const int N, 
  !!     const int P, int * K, int * L, float* A, const int LDA, float* B, const int LDB, 
  !!     float* ALPHA, float* BETA, float* U, const int LDU, float* V, const int LDV, float* Q, 
  !!     const int LDQ, float* WORK, int* IWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GGSVD_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBU, 
  !!     const char JOBV, const char JOBQ, const int M, const int N, const int P, int * K, 
  !!     int * L, float * A, const int LDA, float * B, const int LDB, float * ALPHA, float * BETA, 
  !!     float * U, const int LDU, float * V, const int LDV, float * Q, const int LDQ, 
  !!     float * WORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GGSVD_float ( selfID, JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, &
        LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, IWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GGSVD_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBU
    character(kind=c_char)      ,intent(in)   ,value              :: JOBV
    character(kind=c_char)      ,intent(in)   ,value              :: JOBQ
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: P
    integer(c_int)                                  ,dimension(*) :: K
    integer(c_int)                                  ,dimension(*) :: L
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)                                   ,dimension(*) :: ALPHA
    real(c_float)                                   ,dimension(*) :: BETA
    real(c_float)                                   ,dimension(*) :: U
    integer(c_int)              ,intent(in)   ,value              :: LDU
    real(c_float)                                   ,dimension(*) :: V
    integer(c_int)              ,intent(in)   ,value              :: LDV
    real(c_float)                                   ,dimension(*) :: Q
    integer(c_int)              ,intent(in)   ,value              :: LDQ
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEEV(const char JOBVL, const char JOBVR, const int N, double* A, const int LDA, 
  !!     double* WR, double* WI, double* VL, const int LDVL, double* VR, const int LDVR, 
  !!     double* WORK, const int LWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEEV_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBVL, 
  !!     const char JOBVR, const int N, double * A, const int LDA, double * WR, double * WI, 
  !!     double * VL, const int LDVL, double * VR, const int LDVR, double * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_GEEV_double ( selfID, JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, &
        VR, LDVR, WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_GEEV_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVL
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVR
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: WR
    real(c_double)                                  ,dimension(*) :: WI
    real(c_double)                                  ,dimension(*) :: VL
    integer(c_int)              ,intent(in)   ,value              :: LDVL
    real(c_double)                                  ,dimension(*) :: VR
    integer(c_int)              ,intent(in)   ,value              :: LDVR
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEEV(const char JOBVL, const char JOBVR, const int N, float* A, const int LDA, float* WR, 
  !!     float* WI, float* VL, const int LDVL, float* VR, const int LDVR, float* WORK, 
  !!     const int LWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEEV_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBVL, 
  !!     const char JOBVR, const int N, float * A, const int LDA, float * WR, float * WI, 
  !!     float * VL, const int LDVL, float * VR, const int LDVR, float * WORK, const int LWORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_GEEV_float ( selfID, JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, &
        VR, LDVR, WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_GEEV_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVL
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVR
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: WR
    real(c_float)                                   ,dimension(*) :: WI
    real(c_float)                                   ,dimension(*) :: VL
    integer(c_int)              ,intent(in)   ,value              :: LDVL
    real(c_float)                                   ,dimension(*) :: VR
    integer(c_int)              ,intent(in)   ,value              :: LDVR
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SPEV(const char JOBZ, const char UPLO, const int N, double* AP, double* W, double* Z, 
  !!     int LDZ, double* WORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SPEV_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char UPLO, const int N, double * AP, double * W, double * Z, int LDZ, 
  !!     double * WORK, int * INFO );

  subroutine Epetra_LAPACK_SPEV_double ( selfID, JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO ) &
        bind(C,name='Epetra_LAPACK_SPEV_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: AP
    real(c_double)                                  ,dimension(*) :: W
    real(c_double)                                  ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SPEV(const char JOBZ, const char UPLO, const int N, float* AP, float* W, float* Z, 
  !!     int LDZ, float* WORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SPEV_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char UPLO, const int N, float * AP, float * W, float * Z, int LDZ, float * WORK, 
  !!     int * INFO );

  subroutine Epetra_LAPACK_SPEV_float ( selfID, JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO ) &
        bind(C,name='Epetra_LAPACK_SPEV_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: AP
    real(c_float)                                   ,dimension(*) :: W
    real(c_float)                                   ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SPGV(const int ITYPE, const char JOBZ, const char UPLO, const int N, double* AP, 
  !!     double* BP, double* W, double* Z, const int LDZ, double* WORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SPGV_double ( CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, 
  !!     const char JOBZ, const char UPLO, const int N, double * AP, double * BP, double * W, 
  !!     double * Z, const int LDZ, double * WORK, int * INFO );

  subroutine Epetra_LAPACK_SPGV_double ( selfID, ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, &
        WORK, INFO ) bind(C,name='Epetra_LAPACK_SPGV_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_char ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: ITYPE
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: AP
    real(c_double)                                  ,dimension(*) :: BP
    real(c_double)                                  ,dimension(*) :: W
    real(c_double)                                  ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SPGV(const int ITYPE, const char JOBZ, const char UPLO, const int N, float* AP, 
  !!     float* BP, float* W, float* Z, const int LDZ, float* WORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SPGV_float ( CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, 
  !!     const char JOBZ, const char UPLO, const int N, float * AP, float * BP, float * W, 
  !!     float * Z, const int LDZ, float * WORK, int * INFO );

  subroutine Epetra_LAPACK_SPGV_float ( selfID, ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, &
        WORK, INFO ) bind(C,name='Epetra_LAPACK_SPGV_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_char ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: ITYPE
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: AP
    real(c_float)                                   ,dimension(*) :: BP
    real(c_float)                                   ,dimension(*) :: W
    real(c_float)                                   ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYEV(const char JOBZ, const char UPLO, const int N, double* A, const int LDA, double* W, 
  !!     double* WORK, const int LWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYEV_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char UPLO, const int N, double * A, const int LDA, double * W, double * WORK, 
  !!     const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_SYEV_double ( selfID, JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, &
        INFO ) bind(C,name='Epetra_LAPACK_SYEV_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: W
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYEV(const char JOBZ, const char UPLO, const int N, float* A, const int LDA, float* W, 
  !!     float* WORK, const int LWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYEV_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char UPLO, const int N, float * A, const int LDA, float * W, float * WORK, 
  !!     const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_SYEV_float ( selfID, JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_SYEV_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: W
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYEVD(const char JOBZ, const char UPLO, const int N, double* A, const int LDA, double* W, 
  !!     double* WORK, const int LWORK, int* IWORK, const int LIWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYEVD_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char UPLO, const int N, double * A, const int LDA, double * W, double * WORK, 
  !!     const int LWORK, int * IWORK, const int LIWORK, int * INFO );

  subroutine Epetra_LAPACK_SYEVD_double ( selfID, JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, &
        IWORK, LIWORK, INFO ) bind(C,name='Epetra_LAPACK_SYEVD_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: W
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)              ,intent(in)   ,value              :: LIWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYEVD(const char JOBZ, const char UPLO, const int N, float* A, const int LDA, float* W, 
  !!     float* WORK, const int LWORK, int* IWORK, const int LIWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYEVD_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char UPLO, const int N, float * A, const int LDA, float * W, float * WORK, 
  !!     const int LWORK, int * IWORK, const int LIWORK, int * INFO );

  subroutine Epetra_LAPACK_SYEVD_float ( selfID, JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, &
        IWORK, LIWORK, INFO ) bind(C,name='Epetra_LAPACK_SYEVD_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: W
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)              ,intent(in)   ,value              :: LIWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYEVX(const char JOBZ, const char RANGE, const char UPLO, const int N, double* A, 
  !!     const int LDA, const double* VL, const double* VU, const int* IL, const int* IU, 
  !!     const double ABSTOL, int * M, double* W, double* Z, const int LDZ, double* WORK, 
  !!     const int LWORK, int* IWORK, int* IFAIL, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYEVX_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char RANGE, const char UPLO, const int N, double * A, const int LDA, 
  !!     const double * VL, const double * VU, const int * IL, const int * IU, 
  !!     const double ABSTOL, int * M, double * W, double * Z, const int LDZ, double * WORK, 
  !!     const int LWORK, int * IWORK, int * IFAIL, int * INFO );

  subroutine Epetra_LAPACK_SYEVX_double ( selfID, JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, &
        IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO ) &
        bind(C,name='Epetra_LAPACK_SYEVX_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: RANGE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)         ,dimension(*) :: VL
    real(c_double)              ,intent(in)         ,dimension(*) :: VU
    integer(c_int)              ,intent(in)         ,dimension(*) :: IL
    integer(c_int)              ,intent(in)         ,dimension(*) :: IU
    real(c_double)              ,intent(in)   ,value              :: ABSTOL
    integer(c_int)                                  ,dimension(*) :: M
    real(c_double)                                  ,dimension(*) :: W
    real(c_double)                                  ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: IFAIL
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYEVX(const char JOBZ, const char RANGE, const char UPLO, const int N, float* A, 
  !!     const int LDA, const float* VL, const float* VU, const int* IL, const int* IU, 
  !!     const float ABSTOL, int * M, float* W, float* Z, const int LDZ, float* WORK, 
  !!     const int LWORK, int* IWORK, int* IFAIL, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYEVX_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char RANGE, const char UPLO, const int N, float * A, const int LDA, 
  !!     const float * VL, const float * VU, const int * IL, const int * IU, const float ABSTOL, 
  !!     int * M, float * W, float * Z, const int LDZ, float * WORK, const int LWORK, int * IWORK, 
  !!     int * IFAIL, int * INFO );

  subroutine Epetra_LAPACK_SYEVX_float ( selfID, JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, &
        IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO ) &
        bind(C,name='Epetra_LAPACK_SYEVX_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: RANGE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)         ,dimension(*) :: VL
    real(c_float)               ,intent(in)         ,dimension(*) :: VU
    integer(c_int)              ,intent(in)         ,dimension(*) :: IL
    integer(c_int)              ,intent(in)         ,dimension(*) :: IU
    real(c_float)               ,intent(in)   ,value              :: ABSTOL
    integer(c_int)                                  ,dimension(*) :: M
    real(c_float)                                   ,dimension(*) :: W
    real(c_float)                                   ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: IFAIL
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYGV(const int ITYPE, const char JOBZ, const char UPLO, const int N, double* A, 
  !!     const int LDA, double* B, const int LDB, double* W, double* WORK, const int LWORK, 
  !!     int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYGV_double ( CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, 
  !!     const char JOBZ, const char UPLO, const int N, double * A, const int LDA, double * B, 
  !!     const int LDB, double * W, double * WORK, const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_SYGV_double ( selfID, ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, &
        WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_SYGV_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_char ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: ITYPE
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)                                  ,dimension(*) :: W
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYGV(const int ITYPE, const char JOBZ, const char UPLO, const int N, float* A, 
  !!     const int LDA, float* B, const int LDB, float* W, float* WORK, const int LWORK, 
  !!     int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYGV_float ( CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, 
  !!     const char JOBZ, const char UPLO, const int N, float * A, const int LDA, float * B, 
  !!     const int LDB, float * W, float * WORK, const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_SYGV_float ( selfID, ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, &
        WORK, LWORK, INFO ) bind(C,name='Epetra_LAPACK_SYGV_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_char ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: ITYPE
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)                                   ,dimension(*) :: W
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYGVX(const int ITYPE, const char JOBZ, const char RANGE, const char UPLO, const int N, 
  !!     double* A, const int LDA, double* B, const int LDB, const double* VL, const double* VU, 
  !!     const int* IL, const int* IU, const double ABSTOL, int* M, double* W, double* Z, 
  !!     const int LDZ, double* WORK, const int LWORK, int* IWORK, int* IFAIL, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYGVX_double ( CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, 
  !!     const char JOBZ, const char RANGE, const char UPLO, const int N, double * A, 
  !!     const int LDA, double * B, const int LDB, const double * VL, const double * VU, 
  !!     const int * IL, const int * IU, const double ABSTOL, int * M, double * W, double * Z, 
  !!     const int LDZ, double * WORK, const int LWORK, int * IWORK, int * IFAIL, int * INFO );

  subroutine Epetra_LAPACK_SYGVX_double ( selfID, ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, &
        LDB, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO ) &
        bind(C,name='Epetra_LAPACK_SYGVX_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_char ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: ITYPE
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: RANGE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)              ,intent(in)         ,dimension(*) :: VL
    real(c_double)              ,intent(in)         ,dimension(*) :: VU
    integer(c_int)              ,intent(in)         ,dimension(*) :: IL
    integer(c_int)              ,intent(in)         ,dimension(*) :: IU
    real(c_double)              ,intent(in)   ,value              :: ABSTOL
    integer(c_int)                                  ,dimension(*) :: M
    real(c_double)                                  ,dimension(*) :: W
    real(c_double)                                  ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: IFAIL
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYGVX(const int ITYPE, const char JOBZ, const char RANGE, const char UPLO, const int N, 
  !!     float* A, const int LDA, float* B, const int LDB, const float* VL, const float* VU, 
  !!     const int* IL, const int* IU, const float ABSTOL, int* M, float* W, float* Z, 
  !!     const int LDZ, float* WORK, const int LWORK, int* IWORK, int* IFAIL, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYGVX_float ( CT_Epetra_LAPACK_ID_t selfID, const int ITYPE, 
  !!     const char JOBZ, const char RANGE, const char UPLO, const int N, float * A, 
  !!     const int LDA, float * B, const int LDB, const float * VL, const float * VU, 
  !!     const int * IL, const int * IU, const float ABSTOL, int * M, float * W, float * Z, 
  !!     const int LDZ, float * WORK, const int LWORK, int * IWORK, int * IFAIL, int * INFO );

  subroutine Epetra_LAPACK_SYGVX_float ( selfID, ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, &
        LDB, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO ) &
        bind(C,name='Epetra_LAPACK_SYGVX_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_char ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: ITYPE
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: RANGE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)               ,intent(in)         ,dimension(*) :: VL
    real(c_float)               ,intent(in)         ,dimension(*) :: VU
    integer(c_int)              ,intent(in)         ,dimension(*) :: IL
    integer(c_int)              ,intent(in)         ,dimension(*) :: IU
    real(c_float)               ,intent(in)   ,value              :: ABSTOL
    integer(c_int)                                  ,dimension(*) :: M
    real(c_float)                                   ,dimension(*) :: W
    real(c_float)                                   ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: IFAIL
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYEVR(const char JOBZ, const char RANGE, const char UPLO, const int N, double* A, 
  !!     const int LDA, const double* VL, const double* VU, const int *IL, const int *IU, 
  !!     const double ABSTOL, int* M, double* W, double* Z, const int LDZ, int* ISUPPZ, 
  !!     double* WORK, const int LWORK, int* IWORK, const int LIWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYEVR_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char RANGE, const char UPLO, const int N, double * A, const int LDA, 
  !!     const double * VL, const double * VU, const int * IL, const int * IU, 
  !!     const double ABSTOL, int * M, double * W, double * Z, const int LDZ, int * ISUPPZ, 
  !!     double * WORK, const int LWORK, int * IWORK, const int LIWORK, int * INFO );

  subroutine Epetra_LAPACK_SYEVR_double ( selfID, JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, &
        IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_SYEVR_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: RANGE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)              ,intent(in)         ,dimension(*) :: VL
    real(c_double)              ,intent(in)         ,dimension(*) :: VU
    integer(c_int)              ,intent(in)         ,dimension(*) :: IL
    integer(c_int)              ,intent(in)         ,dimension(*) :: IU
    real(c_double)              ,intent(in)   ,value              :: ABSTOL
    integer(c_int)                                  ,dimension(*) :: M
    real(c_double)                                  ,dimension(*) :: W
    real(c_double)                                  ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    integer(c_int)                                  ,dimension(*) :: ISUPPZ
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)              ,intent(in)   ,value              :: LIWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void SYEVR(const char JOBZ, const char RANGE, const char UPLO, const int N, float* A, 
  !!     const int LDA, const float* VL, const float* VU, const int *IL, const int *IU, 
  !!     const float ABSTOL, int* M, float* W, float* Z, const int LDZ, int* ISUPPZ, float* WORK, 
  !!     const int LWORK, int* IWORK, const int LIWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_SYEVR_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, 
  !!     const char RANGE, const char UPLO, const int N, float * A, const int LDA, 
  !!     const float * VL, const float * VU, const int * IL, const int * IU, const float ABSTOL, 
  !!     int * M, float * W, float * Z, const int LDZ, int * ISUPPZ, float * WORK, 
  !!     const int LWORK, int * IWORK, const int LIWORK, int * INFO );

  subroutine Epetra_LAPACK_SYEVR_float ( selfID, JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, &
        IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_SYEVR_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    character(kind=c_char)      ,intent(in)   ,value              :: RANGE
    character(kind=c_char)      ,intent(in)   ,value              :: UPLO
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)               ,intent(in)         ,dimension(*) :: VL
    real(c_float)               ,intent(in)         ,dimension(*) :: VU
    integer(c_int)              ,intent(in)         ,dimension(*) :: IL
    integer(c_int)              ,intent(in)         ,dimension(*) :: IU
    real(c_float)               ,intent(in)   ,value              :: ABSTOL
    integer(c_int)                                  ,dimension(*) :: M
    real(c_float)                                   ,dimension(*) :: W
    real(c_float)                                   ,dimension(*) :: Z
    integer(c_int)              ,intent(in)   ,value              :: LDZ
    integer(c_int)                                  ,dimension(*) :: ISUPPZ
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)              ,intent(in)   ,value              :: LIWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, 
  !!     const int N, double* A, const int LDA, double* WR, double* WI, double* VL, 
  !!     const int LDVL, double* VR, const int LDVR, int* ILO, int* IHI, double* SCALE, 
  !!     double* ABNRM, double* RCONDE, double* RCONDV, double* WORK, const int LWORK, int* IWORK, 
  !!     int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEEVX_double ( CT_Epetra_LAPACK_ID_t selfID, const char BALANC, 
  !!     const char JOBVL, const char JOBVR, const char SENSE, const int N, double * A, 
  !!     const int LDA, double * WR, double * WI, double * VL, const int LDVL, double * VR, 
  !!     const int LDVR, int * ILO, int * IHI, double * SCALE, double * ABNRM, double * RCONDE, 
  !!     double * RCONDV, double * WORK, const int LWORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GEEVX_double ( selfID, BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, &
        WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, &
        IWORK, INFO ) bind(C,name='Epetra_LAPACK_GEEVX_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: BALANC
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVL
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVR
    character(kind=c_char)      ,intent(in)   ,value              :: SENSE
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: WR
    real(c_double)                                  ,dimension(*) :: WI
    real(c_double)                                  ,dimension(*) :: VL
    integer(c_int)              ,intent(in)   ,value              :: LDVL
    real(c_double)                                  ,dimension(*) :: VR
    integer(c_int)              ,intent(in)   ,value              :: LDVR
    integer(c_int)                                  ,dimension(*) :: ILO
    integer(c_int)                                  ,dimension(*) :: IHI
    real(c_double)                                  ,dimension(*) :: SCALE
    real(c_double)                                  ,dimension(*) :: ABNRM
    real(c_double)                                  ,dimension(*) :: RCONDE
    real(c_double)                                  ,dimension(*) :: RCONDV
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, 
  !!     const int N, float* A, const int LDA, float* WR, float* WI, float* VL, const int LDVL, 
  !!     float* VR, const int LDVR, int* ILO, int* IHI, float* SCALE, float* ABNRM, float* RCONDE, 
  !!     float* RCONDV, float* WORK, const int LWORK, int* IWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GEEVX_float ( CT_Epetra_LAPACK_ID_t selfID, const char BALANC, 
  !!     const char JOBVL, const char JOBVR, const char SENSE, const int N, float * A, 
  !!     const int LDA, float * WR, float * WI, float * VL, const int LDVL, float * VR, 
  !!     const int LDVR, int * ILO, int * IHI, float * SCALE, float * ABNRM, float * RCONDE, 
  !!     float * RCONDV, float * WORK, const int LWORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GEEVX_float ( selfID, BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, &
        WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, &
        INFO ) bind(C,name='Epetra_LAPACK_GEEVX_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: BALANC
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVL
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVR
    character(kind=c_char)      ,intent(in)   ,value              :: SENSE
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: WR
    real(c_float)                                   ,dimension(*) :: WI
    real(c_float)                                   ,dimension(*) :: VL
    integer(c_int)              ,intent(in)   ,value              :: LDVL
    real(c_float)                                   ,dimension(*) :: VR
    integer(c_int)              ,intent(in)   ,value              :: LDVR
    integer(c_int)                                  ,dimension(*) :: ILO
    integer(c_int)                                  ,dimension(*) :: IHI
    real(c_float)                                   ,dimension(*) :: SCALE
    real(c_float)                                   ,dimension(*) :: ABNRM
    real(c_float)                                   ,dimension(*) :: RCONDE
    real(c_float)                                   ,dimension(*) :: RCONDV
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GESDD(const char JOBZ, const int M, const int N, double* A, const int LDA, double* S, 
  !!     double* U, const int LDU, double* VT, const int LDVT, double* WORK, const int LWORK, 
  !!     int* IWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GESDD_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const int M, 
  !!     const int N, double * A, const int LDA, double * S, double * U, const int LDU, 
  !!     double * VT, const int LDVT, double * WORK, const int LWORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GESDD_double ( selfID, JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, &
        WORK, LWORK, IWORK, INFO ) bind(C,name='Epetra_LAPACK_GESDD_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: S
    real(c_double)                                  ,dimension(*) :: U
    integer(c_int)              ,intent(in)   ,value              :: LDU
    real(c_double)                                  ,dimension(*) :: VT
    integer(c_int)              ,intent(in)   ,value              :: LDVT
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GESDD(const char JOBZ, const int M, const int N, float* A, const int LDA, float* S, 
  !!     float* U, const int LDU, float* VT, const int LDVT, float* WORK, const int LWORK, 
  !!     int* IWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GESDD_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBZ, const int M, 
  !!     const int N, float * A, const int LDA, float * S, float * U, const int LDU, float * VT, 
  !!     const int LDVT, float * WORK, const int LWORK, int * IWORK, int * INFO );

  subroutine Epetra_LAPACK_GESDD_float ( selfID, JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, &
        WORK, LWORK, IWORK, INFO ) bind(C,name='Epetra_LAPACK_GESDD_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBZ
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: S
    real(c_float)                                   ,dimension(*) :: U
    integer(c_int)              ,intent(in)   ,value              :: LDU
    real(c_float)                                   ,dimension(*) :: VT
    integer(c_int)              ,intent(in)   ,value              :: LDVT
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: IWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GGEV(const char JOBVL, const char JOBVR, const int N, double* A, const int LDA, 
  !!     double* B, const int LDB, double* ALPHAR, double* ALPHAI, double* BETA, double* VL, 
  !!     const int LDVL, double* VR, const int LDVR, double* WORK, const int LWORK, 
  !!     int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GGEV_double ( CT_Epetra_LAPACK_ID_t selfID, const char JOBVL, 
  !!     const char JOBVR, const int N, double * A, const int LDA, double * B, const int LDB, 
  !!     double * ALPHAR, double * ALPHAI, double * BETA, double * VL, const int LDVL, 
  !!     double * VR, const int LDVR, double * WORK, const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_GGEV_double ( selfID, JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, &
        ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GGEV_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVL
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVR
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)                                  ,dimension(*) :: ALPHAR
    real(c_double)                                  ,dimension(*) :: ALPHAI
    real(c_double)                                  ,dimension(*) :: BETA
    real(c_double)                                  ,dimension(*) :: VL
    integer(c_int)              ,intent(in)   ,value              :: LDVL
    real(c_double)                                  ,dimension(*) :: VR
    integer(c_int)              ,intent(in)   ,value              :: LDVR
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GGEV(const char JOBVL, const char JOBVR, const int N, float* A, const int LDA, float* B, 
  !!     const int LDB, float* ALPHAR, float* ALPHAI, float* BETA, float* VL, const int LDVL, 
  !!     float* VR, const int LDVR, float* WORK, const int LWORK, int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GGEV_float ( CT_Epetra_LAPACK_ID_t selfID, const char JOBVL, 
  !!     const char JOBVR, const int N, float * A, const int LDA, float * B, const int LDB, 
  !!     float * ALPHAR, float * ALPHAI, float * BETA, float * VL, const int LDVL, float * VR, 
  !!     const int LDVR, float * WORK, const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_GGEV_float ( selfID, JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, &
        ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO ) &
        bind(C,name='Epetra_LAPACK_GGEV_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVL
    character(kind=c_char)      ,intent(in)   ,value              :: JOBVR
    integer(c_int)              ,intent(in)   ,value              :: N
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)                                   ,dimension(*) :: ALPHAR
    real(c_float)                                   ,dimension(*) :: ALPHAI
    real(c_float)                                   ,dimension(*) :: BETA
    real(c_float)                                   ,dimension(*) :: VL
    integer(c_int)              ,intent(in)   ,value              :: LDVL
    real(c_float)                                   ,dimension(*) :: VR
    integer(c_int)              ,intent(in)   ,value              :: LDVR
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GGLSE(const int M, const int N, const int P, double* A, const int LDA, double* B, 
  !!     const int LDB, double* C, double* D, double* X, double* WORK, const int LWORK, 
  !!     int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GGLSE_double ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     const int P, double * A, const int LDA, double * B, const int LDB, double * C, 
  !!     double * D, double * X, double * WORK, const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_GGLSE_double ( selfID, M, N, P, A, LDA, B, LDB, C, D, X, WORK, &
        LWORK, INFO ) bind(C,name='Epetra_LAPACK_GGLSE_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: P
    real(c_double)                                  ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_double)                                  ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_double)                                  ,dimension(*) :: C
    real(c_double)                                  ,dimension(*) :: D
    real(c_double)                                  ,dimension(*) :: X
    real(c_double)                                  ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void GGLSE(const int M, const int N, const int P, float* A, const int LDA, float* B, 
  !!     const int LDB, float* C, float* D, float* X, float* WORK, const int LWORK, 
  !!     int* INFO) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_GGLSE_float ( CT_Epetra_LAPACK_ID_t selfID, const int M, const int N, 
  !!     const int P, float * A, const int LDA, float * B, const int LDB, float * C, float * D, 
  !!     float * X, float * WORK, const int LWORK, int * INFO );

  subroutine Epetra_LAPACK_GGLSE_float ( selfID, M, N, P, A, LDA, B, LDB, C, D, X, WORK, &
        LWORK, INFO ) bind(C,name='Epetra_LAPACK_GGLSE_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_int ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: M
    integer(c_int)              ,intent(in)   ,value              :: N
    integer(c_int)              ,intent(in)   ,value              :: P
    real(c_float)                                   ,dimension(*) :: A
    integer(c_int)              ,intent(in)   ,value              :: LDA
    real(c_float)                                   ,dimension(*) :: B
    integer(c_int)              ,intent(in)   ,value              :: LDB
    real(c_float)                                   ,dimension(*) :: C
    real(c_float)                                   ,dimension(*) :: D
    real(c_float)                                   ,dimension(*) :: X
    real(c_float)                                   ,dimension(*) :: WORK
    integer(c_int)              ,intent(in)   ,value              :: LWORK
    integer(c_int)                                  ,dimension(*) :: INFO
  end subroutine


  !> <BR> Original C++ prototype:
  !! void LAMCH ( const char CMACH, float & T) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_LAMCH_float ( CT_Epetra_LAPACK_ID_t selfID, const char CMACH, float * T );

  subroutine Epetra_LAPACK_LAMCH_float ( selfID, CMACH, T ) &
        bind(C,name='Epetra_LAPACK_LAMCH_float')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_float
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: CMACH
    real(c_float)               ,intent(inout)                    :: T
  end subroutine


  !> <BR> Original C++ prototype:
  !! void LAMCH ( const char CMACH, double & T) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_LAPACK_LAMCH_double ( CT_Epetra_LAPACK_ID_t selfID, const char CMACH, double * T );

  subroutine Epetra_LAPACK_LAMCH_double ( selfID, CMACH, T ) &
        bind(C,name='Epetra_LAPACK_LAMCH_double')
    import :: FT_Epetra_LAPACK_ID_t ,c_char ,c_double
    
    type(FT_Epetra_LAPACK_ID_t) ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)   ,value              :: CMACH
    real(c_double)              ,intent(inout)                    :: T
  end subroutine


!> @}


!> @name Epetra_FECrsMatrix interface
!! @{

  ! _________________ Epetra_FECrsMatrix interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_FECrsMatrix_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_Degeneralize')
    import :: FT_Epetra_FECrsMatrix_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_FECrsMatrix_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_FECrsMatrix_Generalize ( CT_Epetra_FECrsMatrix_ID_t id );

  function Epetra_FECrsMatrix_Generalize ( id ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_FECrsMatrix_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                  :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int* NumEntriesPerRow, 
  !!     bool ignoreNonLocalEntries=false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_Var ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_Map_ID_t RowMapID, int * NumEntriesPerRow, boolean ignoreNonLocalEntries );

  function Epetra_FECrsMatrix_Create_Var ( CV, RowMapID, NumEntriesPerRow, &
        ignoreNonLocalEntries ) result(that) bind(C,name='Epetra_FECrsMatrix_Create_Var')
    import :: FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_Map_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_FECrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: RowMapID
    integer(c_int)                                      ,dimension(*) :: NumEntriesPerRow
    integer(FT_boolean_t)           ,intent(in)   ,value              :: ignoreNonLocalEntries
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int NumEntriesPerRow, 
  !!     bool ignoreNonLocalEntries=false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_Map_ID_t RowMapID, int NumEntriesPerRow, boolean ignoreNonLocalEntries );

  function Epetra_FECrsMatrix_Create ( CV, RowMapID, NumEntriesPerRow, &
        ignoreNonLocalEntries ) result(that) bind(C,name='Epetra_FECrsMatrix_Create')
    import :: FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_Map_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_FECrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: RowMapID
    integer(c_int)                  ,intent(in)   ,value              :: NumEntriesPerRow
    integer(FT_boolean_t)           ,intent(in)   ,value              :: ignoreNonLocalEntries
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, 
  !!     int* NumEntriesPerRow, bool ignoreNonLocalEntries=false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_WithColMap_Var ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_Map_ID_t RowMapID, CT_Epetra_Map_ID_t ColMapID, int * NumEntriesPerRow, 
  !!     boolean ignoreNonLocalEntries );

  function Epetra_FECrsMatrix_Create_WithColMap_Var ( CV, RowMapID, ColMapID, &
        NumEntriesPerRow, ignoreNonLocalEntries ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_Create_WithColMap_Var')
    import :: FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_Map_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_FECrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: RowMapID
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: ColMapID
    integer(c_int)                                      ,dimension(*) :: NumEntriesPerRow
    integer(FT_boolean_t)           ,intent(in)   ,value              :: ignoreNonLocalEntries
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const Epetra_Map& ColMap, 
  !!     int NumEntriesPerRow, bool ignoreNonLocalEntries=false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_WithColMap ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_Map_ID_t RowMapID, CT_Epetra_Map_ID_t ColMapID, int NumEntriesPerRow, 
  !!     boolean ignoreNonLocalEntries );

  function Epetra_FECrsMatrix_Create_WithColMap ( CV, RowMapID, ColMapID, NumEntriesPerRow, &
        ignoreNonLocalEntries ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_Create_WithColMap')
    import :: FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,FT_Epetra_Map_ID_t , &
          c_int ,FT_boolean_t
    
    type(FT_Epetra_FECrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: RowMapID
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: ColMapID
    integer(c_int)                  ,intent(in)   ,value              :: NumEntriesPerRow
    integer(FT_boolean_t)           ,intent(in)   ,value              :: ignoreNonLocalEntries
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_FECrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& Graph, 
  !!     bool ignoreNonLocalEntries=false);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Create_FromGraph ( CT_Epetra_DataAccess_E_t CV, 
  !!     CT_Epetra_CrsGraph_ID_t GraphID, boolean ignoreNonLocalEntries );

  function Epetra_FECrsMatrix_Create_FromGraph ( CV, GraphID, ignoreNonLocalEntries ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_Create_FromGraph')
    import :: FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_DataAccess_E_t , &
          FT_Epetra_CrsGraph_ID_t ,FT_boolean_t
    
    type(FT_Epetra_FECrsMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t),intent(in)   ,value              :: CV
    type(FT_Epetra_CrsGraph_ID_t)   ,intent(in)   ,value              :: GraphID
    integer(FT_boolean_t)           ,intent(in)   ,value              :: ignoreNonLocalEntries
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_FECrsMatrix(const Epetra_FECrsMatrix& src);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_FECrsMatrix_ID_t Epetra_FECrsMatrix_Duplicate ( CT_Epetra_FECrsMatrix_ID_t srcID );

  function Epetra_FECrsMatrix_Duplicate ( srcID ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_Duplicate')
    import :: FT_Epetra_FECrsMatrix_ID_t
    
    type(FT_Epetra_FECrsMatrix_ID_t)                                  :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: srcID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_FECrsMatrix();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_FECrsMatrix_Destroy ( CT_Epetra_FECrsMatrix_ID_t * selfID );

  subroutine Epetra_FECrsMatrix_Destroy ( selfID ) &
        bind(C,name='Epetra_FECrsMatrix_Destroy')
    import :: FT_Epetra_FECrsMatrix_ID_t
    
    type(FT_Epetra_FECrsMatrix_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! Epetra_FECrsMatrix& operator=(const Epetra_FECrsMatrix& src);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_FECrsMatrix_Assign ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     CT_Epetra_FECrsMatrix_ID_t srcID );

  subroutine Epetra_FECrsMatrix_Assign ( selfID, srcID ) &
        bind(C,name='Epetra_FECrsMatrix_Assign')
    import :: FT_Epetra_FECrsMatrix_ID_t
    
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: srcID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_SumIntoGlobalValues ( CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, 
  !!     int NumEntries, double * Values, int * Indices );

  function Epetra_FECrsMatrix_SumIntoGlobalValues ( selfID, GlobalRow, NumEntries, Values, &
        Indices ) result(that) bind(C,name='Epetra_FECrsMatrix_SumIntoGlobalValues')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                  ,intent(in)   ,value              :: NumEntries
    real(c_double)                                      ,dimension(*) :: Values
    integer(c_int)                                      ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_InsertGlobalValues ( CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, 
  !!     int NumEntries, double * Values, int * Indices );

  function Epetra_FECrsMatrix_InsertGlobalValues ( selfID, GlobalRow, NumEntries, Values, &
        Indices ) result(that) bind(C,name='Epetra_FECrsMatrix_InsertGlobalValues')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                  ,intent(in)   ,value              :: NumEntries
    real(c_double)                                      ,dimension(*) :: Values
    integer(c_int)                                      ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_ReplaceGlobalValues ( CT_Epetra_FECrsMatrix_ID_t selfID, int GlobalRow, 
  !!     int NumEntries, double * Values, int * Indices );

  function Epetra_FECrsMatrix_ReplaceGlobalValues ( selfID, GlobalRow, NumEntries, Values, &
        Indices ) result(that) bind(C,name='Epetra_FECrsMatrix_ReplaceGlobalValues')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: GlobalRow
    integer(c_int)                  ,intent(in)   ,value              :: NumEntries
    real(c_double)                                      ,dimension(*) :: Values
    integer(c_int)                                      ,dimension(*) :: Indices
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValues(int numIndices, const int* indices, const double* values, 
  !!     int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable_Square ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numIndices, const int * indices, const double * values, int format );

  function Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable_Square ( selfID, numIndices, &
        indices, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable_Square')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numIndices
    integer(c_int)                  ,intent(in)         ,dimension(*) :: indices
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValues(int numRows, const int* rows, int numCols, const int* cols, 
  !!     const double* values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numRows, const int * rows, int numCols, const int * cols, const double * values, 
  !!     int format );

  function Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable ( selfID, numRows, rows, numCols, &
        cols, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_SumIntoGlobalValues_Ftable')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numRows
    integer(c_int)                  ,intent(in)         ,dimension(*) :: rows
    integer(c_int)                  ,intent(in)   ,value              :: numCols
    integer(c_int)                  ,intent(in)         ,dimension(*) :: cols
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValues(int numIndices, const int* indices, const double* const* values, 
  !!     int format=Epetra_FECrsMatrix::ROW_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable_Square ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numIndices, const int * indices, const double* const * values, int format );

  function Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable_Square ( selfID, numIndices, &
        indices, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable_Square')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numIndices
    integer(c_int)                  ,intent(in)         ,dimension(*) :: indices
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValues(int numRows, const int* rows, int numCols, const int* cols, 
  !!     const double* const* values, int format=Epetra_FECrsMatrix::ROW_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numRows, const int * rows, int numCols, const int * cols, 
  !!     const double* const * values, int format );

  function Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable ( selfID, numRows, rows, numCols, &
        cols, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_SumIntoGlobalValues_Ctable')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numRows
    integer(c_int)                  ,intent(in)         ,dimension(*) :: rows
    integer(c_int)                  ,intent(in)   ,value              :: numCols
    integer(c_int)                  ,intent(in)         ,dimension(*) :: cols
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int InsertGlobalValues(int numIndices, const int* indices, const double* values, 
  !!     int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_InsertGlobalValues_Ftable_Square ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numIndices, const int * indices, const double * values, int format );

  function Epetra_FECrsMatrix_InsertGlobalValues_Ftable_Square ( selfID, numIndices, &
        indices, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_InsertGlobalValues_Ftable_Square')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numIndices
    integer(c_int)                  ,intent(in)         ,dimension(*) :: indices
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int InsertGlobalValues(int numRows, const int* rows, int numCols, const int* cols, 
  !!     const double* values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_InsertGlobalValues_Ftable ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numRows, const int * rows, int numCols, const int * cols, const double * values, 
  !!     int format );

  function Epetra_FECrsMatrix_InsertGlobalValues_Ftable ( selfID, numRows, rows, numCols, &
        cols, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_InsertGlobalValues_Ftable')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numRows
    integer(c_int)                  ,intent(in)         ,dimension(*) :: rows
    integer(c_int)                  ,intent(in)   ,value              :: numCols
    integer(c_int)                  ,intent(in)         ,dimension(*) :: cols
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int InsertGlobalValues(int numIndices, const int* indices, const double* const* values, 
  !!     int format=Epetra_FECrsMatrix::ROW_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_InsertGlobalValues_Ctable_Square ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numIndices, const int * indices, const double* const * values, int format );

  function Epetra_FECrsMatrix_InsertGlobalValues_Ctable_Square ( selfID, numIndices, &
        indices, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_InsertGlobalValues_Ctable_Square')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numIndices
    integer(c_int)                  ,intent(in)         ,dimension(*) :: indices
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int InsertGlobalValues(int numRows, const int* rows, int numCols, const int* cols, 
  !!     const double* const* values, int format=Epetra_FECrsMatrix::ROW_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_InsertGlobalValues_Ctable ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numRows, const int * rows, int numCols, const int * cols, 
  !!     const double* const * values, int format );

  function Epetra_FECrsMatrix_InsertGlobalValues_Ctable ( selfID, numRows, rows, numCols, &
        cols, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_InsertGlobalValues_Ctable')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numRows
    integer(c_int)                  ,intent(in)         ,dimension(*) :: rows
    integer(c_int)                  ,intent(in)   ,value              :: numCols
    integer(c_int)                  ,intent(in)         ,dimension(*) :: cols
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValues(int numIndices, const int* indices, const double* values, 
  !!     int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable_Square ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numIndices, const int * indices, const double * values, int format );

  function Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable_Square ( selfID, numIndices, &
        indices, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable_Square')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numIndices
    integer(c_int)                  ,intent(in)         ,dimension(*) :: indices
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValues(int numRows, const int* rows, int numCols, const int* cols, 
  !!     const double* values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numRows, const int * rows, int numCols, const int * cols, const double * values, 
  !!     int format );

  function Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable ( selfID, numRows, rows, numCols, &
        cols, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_ReplaceGlobalValues_Ftable')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numRows
    integer(c_int)                  ,intent(in)         ,dimension(*) :: rows
    integer(c_int)                  ,intent(in)   ,value              :: numCols
    integer(c_int)                  ,intent(in)         ,dimension(*) :: cols
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValues(int numIndices, const int* indices, const double* const* values, 
  !!     int format=Epetra_FECrsMatrix::ROW_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable_Square ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numIndices, const int * indices, const double* const * values, int format );

  function Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable_Square ( selfID, numIndices, &
        indices, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable_Square')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numIndices
    integer(c_int)                  ,intent(in)         ,dimension(*) :: indices
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValues(int numRows, const int* rows, int numCols, const int* cols, 
  !!     const double* const* values, int format=Epetra_FECrsMatrix::ROW_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     int numRows, const int * rows, int numCols, const int * cols, 
  !!     const double* const * values, int format );

  function Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable ( selfID, numRows, rows, numCols, &
        cols, values, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_ReplaceGlobalValues_Ctable')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,c_double
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                  ,intent(in)   ,value              :: numRows
    integer(c_int)                  ,intent(in)         ,dimension(*) :: rows
    integer(c_int)                  ,intent(in)   ,value              :: numCols
    integer(c_int)                  ,intent(in)         ,dimension(*) :: cols
    real(c_double)                  ,intent(in)         ,dimension(*) :: values
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValues(const Epetra_IntSerialDenseVector& indices, 
  !!     const Epetra_SerialDenseMatrix& values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix_Square ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     CT_Epetra_IntSerialDenseVector_ID_t indicesID, CT_Epetra_SerialDenseMatrix_ID_t valuesID, 
  !!     int format );

  function Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix_Square ( selfID, indicesID, &
        valuesID, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix_Square')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_IntSerialDenseVector_ID_t , &
          FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: indicesID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: valuesID
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int SumIntoGlobalValues(const Epetra_IntSerialDenseVector& rows, 
  !!     const Epetra_IntSerialDenseVector& cols, const Epetra_SerialDenseMatrix& values, 
  !!     int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     CT_Epetra_IntSerialDenseVector_ID_t rowsID, CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format );

  function Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix ( selfID, rowsID, colsID, &
        valuesID, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_SumIntoGlobalValues_SubMatrix')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_IntSerialDenseVector_ID_t , &
          FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: rowsID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: colsID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: valuesID
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int InsertGlobalValues(const Epetra_IntSerialDenseVector& indices, 
  !!     const Epetra_SerialDenseMatrix& values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix_Square ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     CT_Epetra_IntSerialDenseVector_ID_t indicesID, CT_Epetra_SerialDenseMatrix_ID_t valuesID, 
  !!     int format );

  function Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix_Square ( selfID, indicesID, &
        valuesID, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix_Square')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_IntSerialDenseVector_ID_t , &
          FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: indicesID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: valuesID
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int InsertGlobalValues(const Epetra_IntSerialDenseVector& rows, 
  !!     const Epetra_IntSerialDenseVector& cols, const Epetra_SerialDenseMatrix& values, 
  !!     int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     CT_Epetra_IntSerialDenseVector_ID_t rowsID, CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format );

  function Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix ( selfID, rowsID, colsID, &
        valuesID, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_InsertGlobalValues_SubMatrix')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_IntSerialDenseVector_ID_t , &
          FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: rowsID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: colsID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: valuesID
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValues(const Epetra_IntSerialDenseVector& indices, 
  !!     const Epetra_SerialDenseMatrix& values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix_Square ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     CT_Epetra_IntSerialDenseVector_ID_t indicesID, CT_Epetra_SerialDenseMatrix_ID_t valuesID, 
  !!     int format );

  function Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix_Square ( selfID, indicesID, &
        valuesID, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix_Square')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_IntSerialDenseVector_ID_t , &
          FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: indicesID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: valuesID
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int ReplaceGlobalValues(const Epetra_IntSerialDenseVector& rows, 
  !!     const Epetra_IntSerialDenseVector& cols, const Epetra_SerialDenseMatrix& values, 
  !!     int format=Epetra_FECrsMatrix::COLUMN_MAJOR);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     CT_Epetra_IntSerialDenseVector_ID_t rowsID, CT_Epetra_IntSerialDenseVector_ID_t colsID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t valuesID, int format );

  function Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix ( selfID, rowsID, colsID, &
        valuesID, format ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_ReplaceGlobalValues_SubMatrix')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_IntSerialDenseVector_ID_t , &
          FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: rowsID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: colsID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: valuesID
    integer(c_int)                  ,intent(in)   ,value              :: format
  end function


  !> <BR> Original C++ prototype:
  !! int GlobalAssemble(bool callFillComplete=true);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_GlobalAssemble ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     boolean callFillComplete );

  function Epetra_FECrsMatrix_GlobalAssemble ( selfID, callFillComplete ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_GlobalAssemble')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,FT_boolean_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)           ,intent(in)   ,value              :: callFillComplete
  end function


  !> <BR> Original C++ prototype:
  !! int GlobalAssemble(const Epetra_Map& domain_map, const Epetra_Map& range_map, 
  !!     bool callFillComplete=true);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_FECrsMatrix_GlobalAssemble_WithMaps ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     CT_Epetra_Map_ID_t domain_mapID, CT_Epetra_Map_ID_t range_mapID, 
  !!     boolean callFillComplete );

  function Epetra_FECrsMatrix_GlobalAssemble_WithMaps ( selfID, domain_mapID, range_mapID, &
        callFillComplete ) result(that) &
        bind(C,name='Epetra_FECrsMatrix_GlobalAssemble_WithMaps')
    import :: c_int ,FT_Epetra_FECrsMatrix_ID_t ,FT_Epetra_Map_ID_t ,FT_boolean_t
    
    integer(c_int)                                                    :: that
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: domain_mapID
    type(FT_Epetra_Map_ID_t)        ,intent(in)   ,value              :: range_mapID
    integer(FT_boolean_t)           ,intent(in)   ,value              :: callFillComplete
  end function


  !> <BR> Original C++ prototype:
  !! void setIgnoreNonLocalEntries(bool flag);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_FECrsMatrix_setIgnoreNonLocalEntries ( CT_Epetra_FECrsMatrix_ID_t selfID, 
  !!     boolean flag );

  subroutine Epetra_FECrsMatrix_setIgnoreNonLocalEntries ( selfID, flag ) &
        bind(C,name='Epetra_FECrsMatrix_setIgnoreNonLocalEntries')
    import :: FT_Epetra_FECrsMatrix_ID_t ,FT_boolean_t
    
    type(FT_Epetra_FECrsMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)           ,intent(in)   ,value              :: flag
  end subroutine


!> @}


!> @name Epetra_IntSerialDenseVector interface
!! @{

  ! _________________ Epetra_IntSerialDenseVector interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_IntSerialDenseVector_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Degeneralize')
    import :: FT_Epetra_IntSerialDenseVector_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_IntSerialDenseVector_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t)         ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_IntSerialDenseVector_Generalize ( CT_Epetra_IntSerialDenseVector_ID_t id );

  function Epetra_IntSerialDenseVector_Generalize ( id ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_IntSerialDenseVector_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                           :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_IntSerialDenseVector();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create_Empty ( );

  function Epetra_IntSerialDenseVector_Create_Empty (  ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Create_Empty')
    import :: FT_Epetra_IntSerialDenseVector_ID_t
    
    type(FT_Epetra_IntSerialDenseVector_ID_t)                                  :: that
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_IntSerialDenseVector(int Length_in);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create ( int Length_in );

  function Epetra_IntSerialDenseVector_Create ( Length_in ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Create')
    import :: FT_Epetra_IntSerialDenseVector_ID_t ,c_int
    
    type(FT_Epetra_IntSerialDenseVector_ID_t)                                  :: that
    integer(c_int)                           ,intent(in)   ,value              :: Length_in
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_IntSerialDenseVector(Epetra_DataAccess CV_in, int* Values_in, int Length_in);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Create_FromArray ( CT_Epetra_DataAccess_E_t CV_in, 
  !!     int * Values_in, int Length_in );

  function Epetra_IntSerialDenseVector_Create_FromArray ( CV_in, Values_in, Length_in ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Create_FromArray')
    import :: FT_Epetra_IntSerialDenseVector_ID_t ,FT_Epetra_DataAccess_E_t ,c_int
    
    type(FT_Epetra_IntSerialDenseVector_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t)        ,intent(in)   ,value              :: CV_in
    integer(c_int)                                               ,dimension(*) :: Values_in
    integer(c_int)                           ,intent(in)   ,value              :: Length_in
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_IntSerialDenseVector(const Epetra_IntSerialDenseVector& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_IntSerialDenseVector_ID_t Epetra_IntSerialDenseVector_Duplicate ( CT_Epetra_IntSerialDenseVector_ID_t SourceID );

  function Epetra_IntSerialDenseVector_Duplicate ( SourceID ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Duplicate')
    import :: FT_Epetra_IntSerialDenseVector_ID_t
    
    type(FT_Epetra_IntSerialDenseVector_ID_t)                                  :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: SourceID
  end function


  !> <BR> Original C++ prototype:
  !! int Size(int Length_in);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_IntSerialDenseVector_Size ( CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  !!     int Length_in );

  function Epetra_IntSerialDenseVector_Size ( selfID, Length_in ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Size')
    import :: c_int ,FT_Epetra_IntSerialDenseVector_ID_t
    
    integer(c_int)                                                             :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                           ,intent(in)   ,value              :: Length_in
  end function


  !> <BR> Original C++ prototype:
  !! int Resize(int Length_in);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_IntSerialDenseVector_Resize ( CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  !!     int Length_in );

  function Epetra_IntSerialDenseVector_Resize ( selfID, Length_in ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Resize')
    import :: c_int ,FT_Epetra_IntSerialDenseVector_ID_t
    
    integer(c_int)                                                             :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                           ,intent(in)   ,value              :: Length_in
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_IntSerialDenseVector ();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_IntSerialDenseVector_Destroy ( CT_Epetra_IntSerialDenseVector_ID_t * selfID );

  subroutine Epetra_IntSerialDenseVector_Destroy ( selfID ) &
        bind(C,name='Epetra_IntSerialDenseVector_Destroy')
    import :: FT_Epetra_IntSerialDenseVector_ID_t
    
    type(FT_Epetra_IntSerialDenseVector_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int& operator () (int Index);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_IntSerialDenseVector_setElement ( CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  !!     int Index, int * value );

  subroutine Epetra_IntSerialDenseVector_setElement ( selfID, Index, value ) &
        bind(C,name='Epetra_IntSerialDenseVector_setElement')
    import :: FT_Epetra_IntSerialDenseVector_ID_t ,c_int
    
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                           ,intent(in)   ,value              :: Index
    integer(c_int)                           ,intent(inout)                    :: value
  end subroutine


  !> <BR> Original C++ prototype:
  !! const int& operator () (int Index) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_IntSerialDenseVector_getElement ( CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  !!     int Index );

  function Epetra_IntSerialDenseVector_getElement ( selfID, Index ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_getElement')
    import :: c_int ,FT_Epetra_IntSerialDenseVector_ID_t
    
    integer(c_int)                                                             :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                           ,intent(in)   ,value              :: Index
  end function


  !> <BR> Original C++ prototype:
  !! int& operator [] (int Index);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_IntSerialDenseVector_setElement_Bracket ( CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  !!     int Index, int * value );

  subroutine Epetra_IntSerialDenseVector_setElement_Bracket ( selfID, Index, value ) &
        bind(C,name='Epetra_IntSerialDenseVector_setElement_Bracket')
    import :: FT_Epetra_IntSerialDenseVector_ID_t ,c_int
    
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                           ,intent(in)   ,value              :: Index
    integer(c_int)                           ,intent(inout)                    :: value
  end subroutine


  !> <BR> Original C++ prototype:
  !! const int& operator [] (int Index) const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_IntSerialDenseVector_getElement_Bracket ( CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  !!     int Index );

  function Epetra_IntSerialDenseVector_getElement_Bracket ( selfID, Index ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_getElement_Bracket')
    import :: c_int ,FT_Epetra_IntSerialDenseVector_ID_t
    
    integer(c_int)                                                             :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                           ,intent(in)   ,value              :: Index
  end function


  !> <BR> Original C++ prototype:
  !! int Random();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_IntSerialDenseVector_Random ( CT_Epetra_IntSerialDenseVector_ID_t selfID );

  function Epetra_IntSerialDenseVector_Random ( selfID ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Random')
    import :: c_int ,FT_Epetra_IntSerialDenseVector_ID_t
    
    integer(c_int)                                                             :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int Length() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_IntSerialDenseVector_Length ( CT_Epetra_IntSerialDenseVector_ID_t selfID );

  function Epetra_IntSerialDenseVector_Length ( selfID ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Length')
    import :: c_int ,FT_Epetra_IntSerialDenseVector_ID_t
    
    integer(c_int)                                                             :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int* Values();
  !> <BR> <BR> CTrilinos prototype:
  !! int * Epetra_IntSerialDenseVector_Values ( CT_Epetra_IntSerialDenseVector_ID_t selfID );

  function Epetra_IntSerialDenseVector_Values ( selfID ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Values')
    import :: c_ptr ,FT_Epetra_IntSerialDenseVector_ID_t
    
    type(c_ptr)                                                                :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! const int* Values() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const int * Epetra_IntSerialDenseVector_Values_Const ( CT_Epetra_IntSerialDenseVector_ID_t selfID );

  function Epetra_IntSerialDenseVector_Values_Const ( selfID ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_Values_Const')
    import :: c_ptr ,FT_Epetra_IntSerialDenseVector_ID_t
    
    type(c_ptr)                                                                :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_DataAccess CV() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_DataAccess_E_t Epetra_IntSerialDenseVector_CV ( CT_Epetra_IntSerialDenseVector_ID_t selfID );

  function Epetra_IntSerialDenseVector_CV ( selfID ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_CV')
    import :: FT_Epetra_DataAccess_E_t ,FT_Epetra_IntSerialDenseVector_ID_t
    
    integer(FT_Epetra_DataAccess_E_t)                                          :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_IntSerialDenseVector& operator = (const Epetra_IntSerialDenseVector& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_IntSerialDenseVector_Assign ( CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  !!     CT_Epetra_IntSerialDenseVector_ID_t SourceID );

  subroutine Epetra_IntSerialDenseVector_Assign ( selfID, SourceID ) &
        bind(C,name='Epetra_IntSerialDenseVector_Assign')
    import :: FT_Epetra_IntSerialDenseVector_ID_t
    
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: SourceID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int MakeViewOf(const Epetra_IntSerialDenseVector& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_IntSerialDenseVector_MakeViewOf ( CT_Epetra_IntSerialDenseVector_ID_t selfID, 
  !!     CT_Epetra_IntSerialDenseVector_ID_t SourceID );

  function Epetra_IntSerialDenseVector_MakeViewOf ( selfID, SourceID ) result(that) &
        bind(C,name='Epetra_IntSerialDenseVector_MakeViewOf')
    import :: c_int ,FT_Epetra_IntSerialDenseVector_ID_t
    
    integer(c_int)                                                             :: that
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_IntSerialDenseVector_ID_t),intent(in)   ,value              :: SourceID
  end function


!> @}


!> @name Epetra_SerialDenseMatrix interface
!! @{

  ! _________________ Epetra_SerialDenseMatrix interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_SerialDenseMatrix_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Degeneralize')
    import :: FT_Epetra_SerialDenseMatrix_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_SerialDenseMatrix_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t)      ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_SerialDenseMatrix_Generalize ( CT_Epetra_SerialDenseMatrix_ID_t id );

  function Epetra_SerialDenseMatrix_Generalize ( id ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_SerialDenseMatrix_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                        :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseMatrix(bool set_object_label=true);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create_Empty ( boolean set_object_label );

  function Epetra_SerialDenseMatrix_Create_Empty ( set_object_label ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Create_Empty')
    import :: FT_Epetra_SerialDenseMatrix_ID_t ,FT_boolean_t
    
    type(FT_Epetra_SerialDenseMatrix_ID_t)                                  :: that
    integer(FT_boolean_t)                 ,intent(in)   ,value              :: set_object_label
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseMatrix(int NumRows, int NumCols, bool set_object_label=true);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create ( int NumRows, int NumCols, 
  !!     boolean set_object_label );

  function Epetra_SerialDenseMatrix_Create ( NumRows, NumCols, set_object_label ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Create')
    import :: FT_Epetra_SerialDenseMatrix_ID_t ,c_int ,FT_boolean_t
    
    type(FT_Epetra_SerialDenseMatrix_ID_t)                                  :: that
    integer(c_int)                        ,intent(in)   ,value              :: NumRows
    integer(c_int)                        ,intent(in)   ,value              :: NumCols
    integer(FT_boolean_t)                 ,intent(in)   ,value              :: set_object_label
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseMatrix(Epetra_DataAccess CV, double* A_in, int LDA_in, int NumRows, 
  !!     int NumCols, bool set_object_label=true);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Create_FromArray ( CT_Epetra_DataAccess_E_t CV, 
  !!     double * A_in, int LDA_in, int NumRows, int NumCols, boolean set_object_label );

  function Epetra_SerialDenseMatrix_Create_FromArray ( CV, A_in, LDA_in, NumRows, NumCols, &
        set_object_label ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Create_FromArray')
    import :: FT_Epetra_SerialDenseMatrix_ID_t ,FT_Epetra_DataAccess_E_t ,c_double ,c_int , &
          FT_boolean_t
    
    type(FT_Epetra_SerialDenseMatrix_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t)     ,intent(in)   ,value              :: CV
    real(c_double)                                            ,dimension(*) :: A_in
    integer(c_int)                        ,intent(in)   ,value              :: LDA_in
    integer(c_int)                        ,intent(in)   ,value              :: NumRows
    integer(c_int)                        ,intent(in)   ,value              :: NumCols
    integer(FT_boolean_t)                 ,intent(in)   ,value              :: set_object_label
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseMatrix(const Epetra_SerialDenseMatrix& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseMatrix_ID_t Epetra_SerialDenseMatrix_Duplicate ( CT_Epetra_SerialDenseMatrix_ID_t SourceID );

  function Epetra_SerialDenseMatrix_Duplicate ( SourceID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Duplicate')
    import :: FT_Epetra_SerialDenseMatrix_ID_t
    
    type(FT_Epetra_SerialDenseMatrix_ID_t)                                  :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: SourceID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_SerialDenseMatrix ();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialDenseMatrix_Destroy ( CT_Epetra_SerialDenseMatrix_ID_t * selfID );

  subroutine Epetra_SerialDenseMatrix_Destroy ( selfID ) &
        bind(C,name='Epetra_SerialDenseMatrix_Destroy')
    import :: FT_Epetra_SerialDenseMatrix_ID_t
    
    type(FT_Epetra_SerialDenseMatrix_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int Shape(int NumRows, int NumCols);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_Shape ( CT_Epetra_SerialDenseMatrix_ID_t selfID, int NumRows, 
  !!     int NumCols );

  function Epetra_SerialDenseMatrix_Shape ( selfID, NumRows, NumCols ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Shape')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                        ,intent(in)   ,value              :: NumRows
    integer(c_int)                        ,intent(in)   ,value              :: NumCols
  end function


  !> <BR> Original C++ prototype:
  !! int Reshape(int NumRows, int NumCols);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_Reshape ( CT_Epetra_SerialDenseMatrix_ID_t selfID, int NumRows, 
  !!     int NumCols );

  function Epetra_SerialDenseMatrix_Reshape ( selfID, NumRows, NumCols ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Reshape')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                        ,intent(in)   ,value              :: NumRows
    integer(c_int)                        ,intent(in)   ,value              :: NumCols
  end function


  !> <BR> Original C++ prototype:
  !! int Multiply(char TransA, char TransB, double ScalarAB, const Epetra_SerialDenseMatrix& A, 
  !!     const Epetra_SerialDenseMatrix& B, double ScalarThis);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_Multiply_Matrix ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     char TransA, char TransB, double ScalarAB, CT_Epetra_SerialDenseMatrix_ID_t AID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t BID, double ScalarThis );

  function Epetra_SerialDenseMatrix_Multiply_Matrix ( selfID, TransA, TransB, ScalarAB, AID, &
        BID, ScalarThis ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Multiply_Matrix')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t ,c_char ,c_double
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                ,intent(in)   ,value              :: TransA
    character(kind=c_char)                ,intent(in)   ,value              :: TransB
    real(c_double)                        ,intent(in)   ,value              :: ScalarAB
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: AID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: BID
    real(c_double)                        ,intent(in)   ,value              :: ScalarThis
  end function


  !> <BR> Original C++ prototype:
  !! int Multiply(bool transA, const Epetra_SerialDenseMatrix& x, Epetra_SerialDenseMatrix& y);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_Multiply_Vector ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     boolean transA, CT_Epetra_SerialDenseMatrix_ID_t xID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t yID );

  function Epetra_SerialDenseMatrix_Multiply_Vector ( selfID, transA, xID, yID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Multiply_Vector')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t ,FT_boolean_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)                 ,intent(in)   ,value              :: transA
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: xID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: yID
  end function


  !> <BR> Original C++ prototype:
  !! int Scale(double ScalarA);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_Scale ( CT_Epetra_SerialDenseMatrix_ID_t selfID, double ScalarA );

  function Epetra_SerialDenseMatrix_Scale ( selfID, ScalarA ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Scale')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t ,c_double
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                        ,intent(in)   ,value              :: ScalarA
  end function


  !> <BR> Original C++ prototype:
  !! virtual double NormOne() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseMatrix_NormOne ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_NormOne ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_NormOne')
    import :: c_double ,FT_Epetra_SerialDenseMatrix_ID_t
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double NormInf() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseMatrix_NormInf ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_NormInf ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_NormInf')
    import :: c_double ,FT_Epetra_SerialDenseMatrix_ID_t
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseMatrix & operator = (const Epetra_SerialDenseMatrix& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialDenseMatrix_Assign ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t SourceID );

  subroutine Epetra_SerialDenseMatrix_Assign ( selfID, SourceID ) &
        bind(C,name='Epetra_SerialDenseMatrix_Assign')
    import :: FT_Epetra_SerialDenseMatrix_ID_t
    
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: SourceID
  end subroutine


  !> <BR> Original C++ prototype:
  !! bool operator==(const Epetra_SerialDenseMatrix& rhs) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_SerialDenseMatrix_IsEqual ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t rhsID );

  function Epetra_SerialDenseMatrix_IsEqual ( selfID, rhsID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_IsEqual')
    import :: FT_boolean_t ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(FT_boolean_t)                                                   :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: rhsID
  end function


  !> <BR> Original C++ prototype:
  !! bool operator!=(const Epetra_SerialDenseMatrix& rhs) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_SerialDenseMatrix_NotEqual ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t rhsID );

  function Epetra_SerialDenseMatrix_NotEqual ( selfID, rhsID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_NotEqual')
    import :: FT_boolean_t ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(FT_boolean_t)                                                   :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: rhsID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseMatrix & operator += (const Epetra_SerialDenseMatrix& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialDenseMatrix_AddTo ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t SourceID );

  subroutine Epetra_SerialDenseMatrix_AddTo ( selfID, SourceID ) &
        bind(C,name='Epetra_SerialDenseMatrix_AddTo')
    import :: FT_Epetra_SerialDenseMatrix_ID_t
    
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: SourceID
  end subroutine


  !> <BR> Original C++ prototype:
  !! double& operator () (int RowIndex, int ColIndex);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialDenseMatrix_setElement ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     int RowIndex, int ColIndex, double * value );

  subroutine Epetra_SerialDenseMatrix_setElement ( selfID, RowIndex, ColIndex, value ) &
        bind(C,name='Epetra_SerialDenseMatrix_setElement')
    import :: FT_Epetra_SerialDenseMatrix_ID_t ,c_int ,c_double
    
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                        ,intent(in)   ,value              :: RowIndex
    integer(c_int)                        ,intent(in)   ,value              :: ColIndex
    real(c_double)                        ,intent(inout)                    :: value
  end subroutine


  !> <BR> Original C++ prototype:
  !! const double& operator () (int RowIndex, int ColIndex) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseMatrix_getElement ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     int RowIndex, int ColIndex );

  function Epetra_SerialDenseMatrix_getElement ( selfID, RowIndex, ColIndex ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_getElement')
    import :: c_double ,FT_Epetra_SerialDenseMatrix_ID_t ,c_int
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                        ,intent(in)   ,value              :: RowIndex
    integer(c_int)                        ,intent(in)   ,value              :: ColIndex
  end function


  !> <BR> Original C++ prototype:
  !! const double* operator [] (int ColIndex) const;
  !> <BR> <BR> CTrilinos prototype:
  !! const double * Epetra_SerialDenseMatrix_getColumn ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     int ColIndex );

  function Epetra_SerialDenseMatrix_getColumn ( selfID, ColIndex ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_getColumn')
    import :: c_ptr ,FT_Epetra_SerialDenseMatrix_ID_t ,c_int
    
    type(c_ptr)                                                             :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                        ,intent(in)   ,value              :: ColIndex
  end function


  !> <BR> Original C++ prototype:
  !! int Random();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_Random ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_Random ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Random')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int M() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_M ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_M ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_M')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int N() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_N ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_N ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_N')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double* A() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double * Epetra_SerialDenseMatrix_A_Const ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_A_Const ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_A_Const')
    import :: c_ptr ,FT_Epetra_SerialDenseMatrix_ID_t
    
    type(c_ptr)                                                             :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double* A();
  !> <BR> <BR> CTrilinos prototype:
  !! double * Epetra_SerialDenseMatrix_A ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_A ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_A')
    import :: c_ptr ,FT_Epetra_SerialDenseMatrix_ID_t
    
    type(c_ptr)                                                             :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int LDA() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_LDA ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_LDA ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_LDA')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_DataAccess CV() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_DataAccess_E_t Epetra_SerialDenseMatrix_CV ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_CV ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_CV')
    import :: FT_Epetra_DataAccess_E_t ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(FT_Epetra_DataAccess_E_t)                                       :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double OneNorm() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseMatrix_OneNorm ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_OneNorm ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_OneNorm')
    import :: c_double ,FT_Epetra_SerialDenseMatrix_ID_t
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double InfNorm() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseMatrix_InfNorm ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_InfNorm ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_InfNorm')
    import :: c_double ,FT_Epetra_SerialDenseMatrix_ID_t
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int SetUseTranspose(bool UseTranspose_in);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_SetUseTranspose ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     boolean UseTranspose_in );

  function Epetra_SerialDenseMatrix_SetUseTranspose ( selfID, UseTranspose_in ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_SetUseTranspose')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t ,FT_boolean_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)                 ,intent(in)   ,value              :: UseTranspose_in
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Apply(const Epetra_SerialDenseMatrix& X, Epetra_SerialDenseMatrix& Y);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_Apply ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t XID, CT_Epetra_SerialDenseMatrix_ID_t YID );

  function Epetra_SerialDenseMatrix_Apply ( selfID, XID, YID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Apply')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ApplyInverse(const Epetra_SerialDenseMatrix & X, Epetra_SerialDenseMatrix & Y);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_ApplyInverse ( CT_Epetra_SerialDenseMatrix_ID_t selfID, 
  !!     CT_Epetra_SerialDenseMatrix_ID_t XID, CT_Epetra_SerialDenseMatrix_ID_t YID );

  function Epetra_SerialDenseMatrix_ApplyInverse ( selfID, XID, YID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_ApplyInverse')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: XID
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const char * Label() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Epetra_SerialDenseMatrix_Label ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_Label ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_Label')
    import :: c_ptr ,FT_Epetra_SerialDenseMatrix_ID_t
    
    type(c_ptr)                                                             :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool UseTranspose() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_SerialDenseMatrix_UseTranspose ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_UseTranspose ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_UseTranspose')
    import :: FT_boolean_t ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(FT_boolean_t)                                                   :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool HasNormInf() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Epetra_SerialDenseMatrix_HasNormInf ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_HasNormInf ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_HasNormInf')
    import :: FT_boolean_t ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(FT_boolean_t)                                                   :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int RowDim() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_RowDim ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_RowDim ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_RowDim')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ColDim() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseMatrix_ColDim ( CT_Epetra_SerialDenseMatrix_ID_t selfID );

  function Epetra_SerialDenseMatrix_ColDim ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseMatrix_ColDim')
    import :: c_int ,FT_Epetra_SerialDenseMatrix_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseMatrix_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name Epetra_SerialDenseVector interface
!! @{

  ! _________________ Epetra_SerialDenseVector interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Epetra_SerialDenseVector_Degeneralize ( id ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Degeneralize')
    import :: FT_Epetra_SerialDenseVector_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Epetra_SerialDenseVector_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t)      ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Epetra_SerialDenseVector_Generalize ( CT_Epetra_SerialDenseVector_ID_t id );

  function Epetra_SerialDenseVector_Generalize ( id ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Epetra_SerialDenseVector_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                        :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseVector();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Create_Empty ( );

  function Epetra_SerialDenseVector_Create_Empty (  ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Create_Empty')
    import :: FT_Epetra_SerialDenseVector_ID_t
    
    type(FT_Epetra_SerialDenseVector_ID_t)                                  :: that
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseVector(int Length);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Create ( int Length );

  function Epetra_SerialDenseVector_Create ( Length ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Create')
    import :: FT_Epetra_SerialDenseVector_ID_t ,c_int
    
    type(FT_Epetra_SerialDenseVector_ID_t)                                  :: that
    integer(c_int)                        ,intent(in)   ,value              :: Length
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseVector(Epetra_DataAccess CV, double* Values, int Length);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Create_FromArray ( CT_Epetra_DataAccess_E_t CV, 
  !!     double * Values, int Length );

  function Epetra_SerialDenseVector_Create_FromArray ( CV, Values, Length ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Create_FromArray')
    import :: FT_Epetra_SerialDenseVector_ID_t ,FT_Epetra_DataAccess_E_t ,c_double ,c_int
    
    type(FT_Epetra_SerialDenseVector_ID_t)                                  :: that
    integer(FT_Epetra_DataAccess_E_t)     ,intent(in)   ,value              :: CV
    real(c_double)                                            ,dimension(*) :: Values
    integer(c_int)                        ,intent(in)   ,value              :: Length
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseVector(const Epetra_SerialDenseVector& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_SerialDenseVector_ID_t Epetra_SerialDenseVector_Duplicate ( CT_Epetra_SerialDenseVector_ID_t SourceID );

  function Epetra_SerialDenseVector_Duplicate ( SourceID ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Duplicate')
    import :: FT_Epetra_SerialDenseVector_ID_t
    
    type(FT_Epetra_SerialDenseVector_ID_t)                                  :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: SourceID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~Epetra_SerialDenseVector ();
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialDenseVector_Destroy ( CT_Epetra_SerialDenseVector_ID_t * selfID );

  subroutine Epetra_SerialDenseVector_Destroy ( selfID ) &
        bind(C,name='Epetra_SerialDenseVector_Destroy')
    import :: FT_Epetra_SerialDenseVector_ID_t
    
    type(FT_Epetra_SerialDenseVector_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! int Size(int Length_in);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseVector_Size ( CT_Epetra_SerialDenseVector_ID_t selfID, int Length_in );

  function Epetra_SerialDenseVector_Size ( selfID, Length_in ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Size')
    import :: c_int ,FT_Epetra_SerialDenseVector_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                        ,intent(in)   ,value              :: Length_in
  end function


  !> <BR> Original C++ prototype:
  !! int Resize(int Length_in);
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseVector_Resize ( CT_Epetra_SerialDenseVector_ID_t selfID, int Length_in );

  function Epetra_SerialDenseVector_Resize ( selfID, Length_in ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Resize')
    import :: c_int ,FT_Epetra_SerialDenseVector_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                        ,intent(in)   ,value              :: Length_in
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_SerialDenseVector& operator = (const Epetra_SerialDenseVector& Source);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialDenseVector_Assign ( CT_Epetra_SerialDenseVector_ID_t selfID, 
  !!     CT_Epetra_SerialDenseVector_ID_t SourceID );

  subroutine Epetra_SerialDenseVector_Assign ( selfID, SourceID ) &
        bind(C,name='Epetra_SerialDenseVector_Assign')
    import :: FT_Epetra_SerialDenseVector_ID_t
    
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: SourceID
  end subroutine


  !> <BR> Original C++ prototype:
  !! double& operator () (int Index);
  !> <BR> <BR> CTrilinos prototype:
  !! void Epetra_SerialDenseVector_setElement ( CT_Epetra_SerialDenseVector_ID_t selfID, int Index, 
  !!     double * value );

  subroutine Epetra_SerialDenseVector_setElement ( selfID, Index, value ) &
        bind(C,name='Epetra_SerialDenseVector_setElement')
    import :: FT_Epetra_SerialDenseVector_ID_t ,c_int ,c_double
    
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                        ,intent(in)   ,value              :: Index
    real(c_double)                        ,intent(inout)                    :: value
  end subroutine


  !> <BR> Original C++ prototype:
  !! const double& operator () (int Index) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseVector_getElement ( CT_Epetra_SerialDenseVector_ID_t selfID, 
  !!     int Index );

  function Epetra_SerialDenseVector_getElement ( selfID, Index ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_getElement')
    import :: c_double ,FT_Epetra_SerialDenseVector_ID_t ,c_int
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                        ,intent(in)   ,value              :: Index
  end function


  !> <BR> Original C++ prototype:
  !! int Random();
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseVector_Random ( CT_Epetra_SerialDenseVector_ID_t selfID );

  function Epetra_SerialDenseVector_Random ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Random')
    import :: c_int ,FT_Epetra_SerialDenseVector_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double Dot(const Epetra_SerialDenseVector & x) const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseVector_Dot ( CT_Epetra_SerialDenseVector_ID_t selfID, 
  !!     CT_Epetra_SerialDenseVector_ID_t xID );

  function Epetra_SerialDenseVector_Dot ( selfID, xID ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Dot')
    import :: c_double ,FT_Epetra_SerialDenseVector_ID_t
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: xID
  end function


  !> <BR> Original C++ prototype:
  !! double Norm1() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseVector_Norm1 ( CT_Epetra_SerialDenseVector_ID_t selfID );

  function Epetra_SerialDenseVector_Norm1 ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Norm1')
    import :: c_double ,FT_Epetra_SerialDenseVector_ID_t
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double Norm2() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseVector_Norm2 ( CT_Epetra_SerialDenseVector_ID_t selfID );

  function Epetra_SerialDenseVector_Norm2 ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Norm2')
    import :: c_double ,FT_Epetra_SerialDenseVector_ID_t
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double NormInf() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double Epetra_SerialDenseVector_NormInf ( CT_Epetra_SerialDenseVector_ID_t selfID );

  function Epetra_SerialDenseVector_NormInf ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_NormInf')
    import :: c_double ,FT_Epetra_SerialDenseVector_ID_t
    
    real(c_double)                                                          :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! int Length() const;
  !> <BR> <BR> CTrilinos prototype:
  !! int Epetra_SerialDenseVector_Length ( CT_Epetra_SerialDenseVector_ID_t selfID );

  function Epetra_SerialDenseVector_Length ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Length')
    import :: c_int ,FT_Epetra_SerialDenseVector_ID_t
    
    integer(c_int)                                                          :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! double* Values() const;
  !> <BR> <BR> CTrilinos prototype:
  !! double * Epetra_SerialDenseVector_Values ( CT_Epetra_SerialDenseVector_ID_t selfID );

  function Epetra_SerialDenseVector_Values ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_Values')
    import :: c_ptr ,FT_Epetra_SerialDenseVector_ID_t
    
    type(c_ptr)                                                             :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! Epetra_DataAccess CV() const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_DataAccess_E_t Epetra_SerialDenseVector_CV ( CT_Epetra_SerialDenseVector_ID_t selfID );

  function Epetra_SerialDenseVector_CV ( selfID ) result(that) &
        bind(C,name='Epetra_SerialDenseVector_CV')
    import :: FT_Epetra_DataAccess_E_t ,FT_Epetra_SerialDenseVector_ID_t
    
    integer(FT_Epetra_DataAccess_E_t)                                       :: that
    type(FT_Epetra_SerialDenseVector_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


  end interface
end module forepetra

