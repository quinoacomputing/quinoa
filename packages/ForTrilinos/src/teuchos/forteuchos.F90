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
module forteuchos
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  use ForTrilinos_enum_wrappers
  implicit none   ! Prevent implicit typing
#ifdef HAVE_MPI
#include <mpif.h>
#endif

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/teuchos/CTeuchos*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

!> @name CommandLineProcessor interface
!! @{

  ! _________________ CommandLineProcessor interface bodies _________________


  !> <BR> Original C++ prototype:
  !! ~CommandLineProcessor();
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_CommandLineProcessor_Destroy ( CT_Teuchos_CommandLineProcessor_ID_t * selfID );

  subroutine Teuchos_CommandLineProcessor_Destroy ( selfID ) &
        bind(C,name='Teuchos_CommandLineProcessor_Destroy')
    import :: FT_Teuchos_CommandLineProcessor_ID_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! CommandLineProcessor( bool throwExceptions = true ,bool recogniseAllOptions = 
  !!     true ,bool addOutputSetupOptions = false );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_CommandLineProcessor_ID_t Teuchos_CommandLineProcessor_Create ( boolean throwExceptions, 
  !!     boolean recogniseAllOptions, boolean addOutputSetupOptions );

  function Teuchos_CommandLineProcessor_Create ( throwExceptions, recogniseAllOptions, &
        addOutputSetupOptions ) result(that) &
        bind(C,name='Teuchos_CommandLineProcessor_Create')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,FT_boolean_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t)                                  :: that
    integer(FT_boolean_t)                     ,intent(in)   ,value              :: throwExceptions
    integer(FT_boolean_t)                     ,intent(in)   ,value              :: recogniseAllOptions
    integer(FT_boolean_t)                     ,intent(in)   ,value              :: addOutputSetupOptions
  end function


  !> <BR> Original C++ prototype:
  !! void throwExceptions( const bool & throwExceptions );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_CommandLineProcessor_throwExceptions_set ( CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  !!     const boolean throwExceptions );

  subroutine Teuchos_CommandLineProcessor_throwExceptions_set ( selfID, throwExceptions ) &
        bind(C,name='Teuchos_CommandLineProcessor_throwExceptions_set')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,FT_boolean_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)                     ,intent(in)   ,value              :: throwExceptions
  end subroutine


  !> <BR> Original C++ prototype:
  !! bool throwExceptions() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_CommandLineProcessor_throwExceptions_get ( CT_Teuchos_CommandLineProcessor_ID_t selfID );

  function Teuchos_CommandLineProcessor_throwExceptions_get ( selfID ) result(that) &
        bind(C,name='Teuchos_CommandLineProcessor_throwExceptions_get')
    import :: FT_boolean_t ,FT_Teuchos_CommandLineProcessor_ID_t
    
    integer(FT_boolean_t)                                                       :: that
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void recogniseAllOptions( const bool & recogniseAllOptions );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_CommandLineProcessor_recogniseAllOptions_set ( CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  !!     const boolean recogniseAllOptions );

  subroutine Teuchos_CommandLineProcessor_recogniseAllOptions_set ( selfID, &
        recogniseAllOptions ) &
        bind(C,name='Teuchos_CommandLineProcessor_recogniseAllOptions_set')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,FT_boolean_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)                     ,intent(in)   ,value              :: recogniseAllOptions
  end subroutine


  !> <BR> Original C++ prototype:
  !! bool recogniseAllOptions() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_CommandLineProcessor_recogniseAllOptions_get ( CT_Teuchos_CommandLineProcessor_ID_t selfID );

  function Teuchos_CommandLineProcessor_recogniseAllOptions_get ( selfID ) result(that) &
        bind(C,name='Teuchos_CommandLineProcessor_recogniseAllOptions_get')
    import :: FT_boolean_t ,FT_Teuchos_CommandLineProcessor_ID_t
    
    integer(FT_boolean_t)                                                       :: that
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void addOutputSetupOptions( const bool &addOutputSetupOptions );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_CommandLineProcessor_addOutputSetupOptions_set ( CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  !!     const boolean addOutputSetupOptions );

  subroutine Teuchos_CommandLineProcessor_addOutputSetupOptions_set ( selfID, &
        addOutputSetupOptions ) &
        bind(C,name='Teuchos_CommandLineProcessor_addOutputSetupOptions_set')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,FT_boolean_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)                     ,intent(in)   ,value              :: addOutputSetupOptions
  end subroutine


  !> <BR> Original C++ prototype:
  !! bool addOutputSetupOptions() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_CommandLineProcessor_addOutputSetupOptions_get ( CT_Teuchos_CommandLineProcessor_ID_t selfID );

  function Teuchos_CommandLineProcessor_addOutputSetupOptions_get ( selfID ) result(that) &
        bind(C,name='Teuchos_CommandLineProcessor_addOutputSetupOptions_get')
    import :: FT_boolean_t ,FT_Teuchos_CommandLineProcessor_ID_t
    
    integer(FT_boolean_t)                                                       :: that
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void setDocString( const char doc_string[] );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_CommandLineProcessor_setDocString ( CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  !!     const char doc_string[] );

  subroutine Teuchos_CommandLineProcessor_setDocString ( selfID, doc_string ) &
        bind(C,name='Teuchos_CommandLineProcessor_setDocString')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: doc_string
  end subroutine


  !> <BR> Original C++ prototype:
  !! void setOption( 
  !!     const char option_true[] ,const char option_false[] ,bool *option_val ,const char documentation[] = 
  !!     NULL );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_CommandLineProcessor_setOption_bool ( CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  !!     const char option_true[], const char option_false[], boolean * option_val, 
  !!     const char documentation[] );

  subroutine Teuchos_CommandLineProcessor_setOption_bool ( selfID, option_true, &
        option_false, option_val, documentation ) &
        bind(C,name='Teuchos_CommandLineProcessor_setOption_bool')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char ,FT_boolean_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_true
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_false
    integer(FT_boolean_t)                                         ,dimension(*) :: option_val
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: documentation
  end subroutine


  !> <BR> Original C++ prototype:
  !! void setOption( const char option_name[] ,int *option_val ,const char documentation[] = 
  !!     NULL ,const bool required = false );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_CommandLineProcessor_setOption_int ( CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  !!     const char option_name[], int * option_val, const char documentation[], 
  !!     const boolean required );

  subroutine Teuchos_CommandLineProcessor_setOption_int ( selfID, option_name, option_val, &
        documentation, required ) bind(C,name='Teuchos_CommandLineProcessor_setOption_int')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char ,c_int ,FT_boolean_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_name
    integer(c_int)                                                ,dimension(*) :: option_val
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: documentation
    integer(FT_boolean_t)                     ,intent(in)   ,value              :: required
  end subroutine


  !> <BR> Original C++ prototype:
  !! void setOption( const char option_name[] ,double *option_val ,const char documentation[] = 
  !!     NULL ,const bool required = false );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_CommandLineProcessor_setOption_double ( CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  !!     const char option_name[], double * option_val, const char documentation[], 
  !!     const boolean required );

  subroutine Teuchos_CommandLineProcessor_setOption_double ( selfID, option_name, &
        option_val, documentation, required ) &
        bind(C,name='Teuchos_CommandLineProcessor_setOption_double')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char ,c_double ,FT_boolean_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_name
    real(c_double)                                                ,dimension(*) :: option_val
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: documentation
    integer(FT_boolean_t)                     ,intent(in)   ,value              :: required
  end subroutine


  !> <BR> Original C++ prototype:
  !! void setOption( 
  !!     const char option_name[] ,std::string *option_val ,const char documentation[] = 
  !!     NULL ,const bool required = false );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_CommandLineProcessor_setOption_str ( CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  !!     const char option_name[], char * option_val[], const char documentation[], 
  !!     const boolean required );

  subroutine Teuchos_CommandLineProcessor_setOption_str ( selfID, option_name, option_val, &
        documentation, required ) bind(C,name='Teuchos_CommandLineProcessor_setOption_str')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char ,FT_boolean_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_name
    character(kind=c_char)                                        ,dimension(*) :: option_val
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: documentation
    integer(FT_boolean_t)                     ,intent(in)   ,value              :: required
  end subroutine


!> @}


!> @name ParameterList interface
!! @{

  ! _________________ ParameterList interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Degeneralize ( CTrilinos_Universal_ID_t id );

  function Teuchos_ParameterList_Degeneralize ( id ) result(that) &
        bind(C,name='Teuchos_ParameterList_Degeneralize')
    import :: FT_Teuchos_ParameterList_ID_t ,ForTrilinos_Universal_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(ForTrilinos_Universal_ID_t)   ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Teuchos_ParameterList_Generalize ( CT_Teuchos_ParameterList_ID_t id );

  function Teuchos_ParameterList_Generalize ( id ) result(that) &
        bind(C,name='Teuchos_ParameterList_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Teuchos_ParameterList_ID_t
    
    type(ForTrilinos_Universal_ID_t)                                     :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create (  );

  function Teuchos_ParameterList_Create (  ) result(that) &
        bind(C,name='Teuchos_ParameterList_Create')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList(const std::string &name);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_WithName ( const char name[] );

  function Teuchos_ParameterList_Create_WithName ( name ) result(that) &
        bind(C,name='Teuchos_ParameterList_Create_WithName')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList(const ParameterList& source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_FromSource ( CT_Teuchos_ParameterList_ID_t sourceID );

  function Teuchos_ParameterList_Create_FromSource ( sourceID ) result(that) &
        bind(C,name='Teuchos_ParameterList_Create_FromSource')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: sourceID
  end function


  !> <BR> Original C++ prototype:
  !! virtual ~ParameterList();
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_ParameterList_Destroy ( CT_Teuchos_ParameterList_ID_t * selfID );

  subroutine Teuchos_ParameterList_Destroy ( selfID ) &
        bind(C,name='Teuchos_ParameterList_Destroy')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! ParameterList& setName( const std::string &name );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setName ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_setName ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_setName')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList& operator=(const ParameterList& source);
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_ParameterList_Assign ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     CT_Teuchos_ParameterList_ID_t sourceID );

  subroutine Teuchos_ParameterList_Assign ( selfID, sourceID ) &
        bind(C,name='Teuchos_ParameterList_Assign')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: sourceID
  end subroutine


  !> <BR> Original C++ prototype:
  !! ParameterList& setParameters(const ParameterList& source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParameters ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     CT_Teuchos_ParameterList_ID_t sourceID );

  function Teuchos_ParameterList_setParameters ( selfID, sourceID ) result(that) &
        bind(C,name='Teuchos_ParameterList_setParameters')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: sourceID
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList& setParametersNotAlreadySet(const ParameterList& source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParametersNotAlreadySet ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     CT_Teuchos_ParameterList_ID_t sourceID );

  function Teuchos_ParameterList_setParametersNotAlreadySet ( selfID, sourceID ) result(that) &
        bind(C,name='Teuchos_ParameterList_setParametersNotAlreadySet')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: sourceID
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList& disableRecursiveValidation();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_disableRecursiveValidation ( CT_Teuchos_ParameterList_ID_t selfID );

  function Teuchos_ParameterList_disableRecursiveValidation ( selfID ) result(that) &
        bind(C,name='Teuchos_ParameterList_disableRecursiveValidation')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> ParameterList& set( std::string const& name, T const& value, 
  !!     std::string const& docString = "" ,RCP<const ParameterEntryValidator> const& validator = 
  !!     Teuchos::null );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_T ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     char const name[], T value, char const docString[] );

  function Teuchos_ParameterList_set_double ( selfID, name, value, docString ) result(that) &
        bind(C,name='Teuchos_ParameterList_set_double')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char ,c_double
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    real(c_double)                     ,intent(in)   ,value              :: value
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: docString
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> ParameterList& set( std::string const& name, T const& value, 
  !!     std::string const& docString = "" ,RCP<const ParameterEntryValidator> const& validator = 
  !!     Teuchos::null );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_T ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     char const name[], T value, char const docString[] );

  function Teuchos_ParameterList_set_int ( selfID, name, value, docString ) result(that) &
        bind(C,name='Teuchos_ParameterList_set_int')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char ,c_int
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    integer(c_int)                     ,intent(in)   ,value              :: value
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: docString
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList& set( std::string const& name, char value[], std::string const& docString = 
  !!     "" ,RCP<const ParameterEntryValidator> const& validator = Teuchos::null );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_str ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     char const name[], char value[], char const docString[] );

  function Teuchos_ParameterList_set_str ( selfID, name, value, docString ) result(that) &
        bind(C,name='Teuchos_ParameterList_set_str')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    character(kind=c_char)                                 ,dimension(*) :: value
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: docString
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList& set( std::string const& name, ParameterList const& value, 
  !!     std::string const& docString = "" );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     char const name[], CT_Teuchos_ParameterList_ID_t valueID, char const docString[] );

  function Teuchos_ParameterList_set ( selfID, name, valueID, docString ) result(that) &
        bind(C,name='Teuchos_ParameterList_set')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: valueID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: docString
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList& setEntry(const std::string& name, const ParameterEntry& entry);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setEntry ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[], CT_Teuchos_ParameterEntry_ID_t entryID );

  function Teuchos_ParameterList_setEntry ( selfID, name, entryID ) result(that) &
        bind(C,name='Teuchos_ParameterList_setEntry')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: entryID
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> T& get(const std::string& name, T def_value);
  !> <BR> <BR> CTrilinos prototype:
  !! T Teuchos_ParameterList_get_def_T ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  !!     T def_value );

  function Teuchos_ParameterList_get_def_double ( selfID, name, def_value ) result(that) &
        bind(C,name='Teuchos_ParameterList_get_def_double')
    import :: c_double ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    real(c_double)                                                       :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    real(c_double)                     ,intent(in)   ,value              :: def_value
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> T& get(const std::string& name, T def_value);
  !> <BR> <BR> CTrilinos prototype:
  !! T Teuchos_ParameterList_get_def_T ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  !!     T def_value );

  function Teuchos_ParameterList_get_def_int ( selfID, name, def_value ) result(that) &
        bind(C,name='Teuchos_ParameterList_get_def_int')
    import :: c_int ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    integer(c_int)                                                       :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    integer(c_int)                     ,intent(in)   ,value              :: def_value
  end function


  !> <BR> Original C++ prototype:
  !! std::string& get(const std::string& name, char def_value[]);
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Teuchos_ParameterList_get_def_char ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[], char def_value[] );

  function Teuchos_ParameterList_get_def_char ( selfID, name, def_value ) result(that) &
        bind(C,name='Teuchos_ParameterList_get_def_char')
    import :: c_char ,FT_Teuchos_ParameterList_ID_t
    
    character(kind=c_char)                                               :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    character(kind=c_char)                                 ,dimension(*) :: def_value
  end function


  !> <BR> Original C++ prototype:
  !! std::string& get(const std::string& name, const char def_value[]);
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Teuchos_ParameterList_get_def_const_char ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[], const char def_value[] );

  function Teuchos_ParameterList_get_def_const_char ( selfID, name, def_value ) result(that) &
        bind(C,name='Teuchos_ParameterList_get_def_const_char')
    import :: c_char ,FT_Teuchos_ParameterList_ID_t
    
    character(kind=c_char)                                               :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: def_value
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> T& get(const std::string& name);
  !> <BR> <BR> CTrilinos prototype:
  !! T Teuchos_ParameterList_get ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  function Teuchos_ParameterList_get_double ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_get_double')
    import :: c_double ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    real(c_double)                                                       :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> T& get(const std::string& name);
  !> <BR> <BR> CTrilinos prototype:
  !! T Teuchos_ParameterList_get ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  function Teuchos_ParameterList_get_int ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_get_int')
    import :: c_int ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    integer(c_int)                                                       :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> const T& get(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! T Teuchos_ParameterList_get_const ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  function Teuchos_ParameterList_get_const_double ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_get_const_double')
    import :: c_double ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    real(c_double)                                                       :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> const T& get(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! T Teuchos_ParameterList_get_const ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  function Teuchos_ParameterList_get_const_int ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_get_const_int')
    import :: c_int ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    integer(c_int)                                                       :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> inline T* getPtr(const std::string& name);
  !> <BR> <BR> CTrilinos prototype:
  !! T * Teuchos_ParameterList_getPtr ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  function Teuchos_ParameterList_getPtr_double ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_getPtr_double')
    import :: c_ptr ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(c_ptr)                                                          :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> inline T* getPtr(const std::string& name);
  !> <BR> <BR> CTrilinos prototype:
  !! T * Teuchos_ParameterList_getPtr ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  function Teuchos_ParameterList_getPtr_int ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_getPtr_int')
    import :: c_ptr ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(c_ptr)                                                          :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> inline const T* getPtr(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! const T * Teuchos_ParameterList_getPtr_const ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_getPtr_const_double ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_getPtr_const_double')
    import :: c_ptr ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(c_ptr)                                                          :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> inline const T* getPtr(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! const T * Teuchos_ParameterList_getPtr_const ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_getPtr_const_int ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_getPtr_const_int')
    import :: c_ptr ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(c_ptr)                                                          :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! ParameterEntry& getEntry(const std::string& name);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_getEntry ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_getEntry')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterEntry_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! inline const ParameterEntry& getEntry(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry_const ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_getEntry_const ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_getEntry_const')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterEntry_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! inline ParameterEntry* getEntryPtr(const std::string& name);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_getEntryPtr ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_getEntryPtr')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterEntry_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! inline const ParameterEntry* getEntryPtr(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr_const ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_getEntryPtr_const ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_getEntryPtr_const')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterEntry_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! bool remove( std::string const& name, bool throwIfNotExists = true );
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterList_remove ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     char const name[], boolean throwIfNotExists );

  function Teuchos_ParameterList_remove ( selfID, name, throwIfNotExists ) result(that) &
        bind(C,name='Teuchos_ParameterList_remove')
    import :: FT_boolean_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    integer(FT_boolean_t)                                                :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    integer(FT_boolean_t)              ,intent(in)   ,value              :: throwIfNotExists
  end function


  !> <BR> Original C++ prototype:
  !! ParameterList& sublist( const std::string& name, bool mustAlreadyExist = 
  !!     false ,const std::string& docString = "" );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[], boolean mustAlreadyExist, const char docString[] );

  function Teuchos_ParameterList_sublist ( selfID, name, mustAlreadyExist, docString ) result(that) &
        bind(C,name='Teuchos_ParameterList_sublist')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char ,FT_boolean_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    integer(FT_boolean_t)              ,intent(in)   ,value              :: mustAlreadyExist
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: docString
  end function


  !> <BR> Original C++ prototype:
  !! const ParameterList& sublist(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist_existing ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_sublist_existing ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_sublist_existing')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! const std::string& name() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Teuchos_ParameterList_name_it ( CT_Teuchos_ParameterList_ID_t selfID );

  function Teuchos_ParameterList_name_it ( selfID ) result(that) &
        bind(C,name='Teuchos_ParameterList_name_it')
    import :: c_char ,FT_Teuchos_ParameterList_ID_t
    
    character(kind=c_char)                                               :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool isParameter(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterList_isParameter ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_isParameter ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_isParameter')
    import :: FT_boolean_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    integer(FT_boolean_t)                                                :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! bool isSublist(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterList_isSublist ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_isSublist ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_isSublist')
    import :: FT_boolean_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    integer(FT_boolean_t)                                                :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> bool isType(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterList_isType ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_isType_double ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_isType_double')
    import :: FT_boolean_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    integer(FT_boolean_t)                                                :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> bool isType(const std::string& name) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterList_isType ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[] );

  function Teuchos_ParameterList_isType_int ( selfID, name ) result(that) &
        bind(C,name='Teuchos_ParameterList_isType_int')
    import :: FT_boolean_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    integer(FT_boolean_t)                                                :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> bool isType(const std::string& name, T* ptr) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterList_isType_type_T ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[], T * ptr );

  function Teuchos_ParameterList_isType_type_double ( selfID, name, ptr ) result(that) &
        bind(C,name='Teuchos_ParameterList_isType_type_double')
    import :: FT_boolean_t ,FT_Teuchos_ParameterList_ID_t ,c_char ,c_double
    
    integer(FT_boolean_t)                                                :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    real(c_double)                                         ,dimension(*) :: ptr
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> bool isType(const std::string& name, T* ptr) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterList_isType_type_T ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     const char name[], T * ptr );

  function Teuchos_ParameterList_isType_type_int ( selfID, name, ptr ) result(that) &
        bind(C,name='Teuchos_ParameterList_isType_type_int')
    import :: FT_boolean_t ,FT_Teuchos_ParameterList_ID_t ,c_char ,c_int
    
    integer(FT_boolean_t)                                                :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    integer(c_int)                                         ,dimension(*) :: ptr
  end function


  !> <BR> Original C++ prototype:
  !! std::string currentParametersString() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Teuchos_ParameterList_currentParametersString ( CT_Teuchos_ParameterList_ID_t selfID );

  function Teuchos_ParameterList_currentParametersString ( selfID ) result(that) &
        bind(C,name='Teuchos_ParameterList_currentParametersString')
    import :: c_char ,FT_Teuchos_ParameterList_ID_t
    
    character(kind=c_char)                                               :: that
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! void validateParameters( ParameterList const& validParamList, int const depth = 1000, 
  !!     EValidateUsed const validateUsed = VALIDATE_USED_ENABLED, 
  !!     EValidateDefaults const validateDefaults = VALIDATE_DEFAULTS_ENABLED ) const;
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_ParameterList_validateParameters ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     CT_Teuchos_ParameterList_ID_t validParamListID, int const depth, 
  !!     const CT_EValidateUsed_E_t validateUsed, 
  !!     const CT_EValidateDefaults_E_t validateDefaults );

  subroutine Teuchos_ParameterList_validateParameters ( selfID, validParamListID, depth, &
        validateUsed, validateDefaults ) &
        bind(C,name='Teuchos_ParameterList_validateParameters')
    import :: FT_Teuchos_ParameterList_ID_t ,c_int ,FT_EValidateUsed_E_t , &
          FT_EValidateDefaults_E_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: validParamListID
    integer(c_int)                     ,intent(in)   ,value              :: depth
    integer(FT_EValidateUsed_E_t)      ,intent(in)   ,value              :: validateUsed
    integer(FT_EValidateDefaults_E_t)  ,intent(in)   ,value              :: validateDefaults
  end subroutine


  !> <BR> Original C++ prototype:
  !! void validateParametersAndSetDefaults( ParameterList const& validParamList, int const depth = 
  !!     1000 );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_ParameterList_validateParametersAndSetDefaults ( CT_Teuchos_ParameterList_ID_t selfID, 
  !!     CT_Teuchos_ParameterList_ID_t validParamListID, int const depth );

  subroutine Teuchos_ParameterList_validateParametersAndSetDefaults ( selfID, &
        validParamListID, depth ) &
        bind(C,name='Teuchos_ParameterList_validateParametersAndSetDefaults')
    import :: FT_Teuchos_ParameterList_ID_t ,c_int
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: validParamListID
    integer(c_int)                     ,intent(in)   ,value              :: depth
  end subroutine


!> @}


!> @name ParameterEntry interface
!! @{

  ! _________________ ParameterEntry interface bodies _________________


  !> <BR> Original C++ prototype:
  !! ParameterEntry();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Create (  );

  function Teuchos_ParameterEntry_Create (  ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_Create')
    import :: FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t)                                  :: that
  end function


  !> <BR> Original C++ prototype:
  !! ParameterEntry(const ParameterEntry& source);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Duplicate ( CT_Teuchos_ParameterEntry_ID_t sourceID );

  function Teuchos_ParameterEntry_Duplicate ( sourceID ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_Duplicate')
    import :: FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t)                                  :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: sourceID
  end function


  !> <BR> Original C++ prototype:
  !! ~ParameterEntry();
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_ParameterEntry_Destroy ( CT_Teuchos_ParameterEntry_ID_t * selfID );

  subroutine Teuchos_ParameterEntry_Destroy ( selfID ) &
        bind(C,name='Teuchos_ParameterEntry_Destroy')
    import :: FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t)                                  :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! ParameterEntry& operator=(const ParameterEntry& source);
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_ParameterEntry_Assign ( CT_Teuchos_ParameterEntry_ID_t selfID, 
  !!     CT_Teuchos_ParameterEntry_ID_t sourceID );

  subroutine Teuchos_ParameterEntry_Assign ( selfID, sourceID ) &
        bind(C,name='Teuchos_ParameterEntry_Assign')
    import :: FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: sourceID
  end subroutine


  !> <BR> Original C++ prototype:
  !! void setAnyValue( const any &value, bool isDefault = false );
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_ParameterEntry_setAnyValue ( CT_Teuchos_ParameterEntry_ID_t selfID, 
  !!     CT_Teuchos_any_ID_t valueID, boolean isDefault );

  subroutine Teuchos_ParameterEntry_setAnyValue ( selfID, valueID, isDefault ) &
        bind(C,name='Teuchos_ParameterEntry_setAnyValue')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_any_ID_t ,FT_boolean_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_any_ID_t)           ,intent(in)   ,value              :: valueID
    integer(FT_boolean_t)               ,intent(in)   ,value              :: isDefault
  end subroutine


  !> <BR> Original C++ prototype:
  !! void setDocString(const std::string &docString);
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_ParameterEntry_setDocString ( CT_Teuchos_ParameterEntry_ID_t selfID, 
  !!     const char docString[] );

  subroutine Teuchos_ParameterEntry_setDocString ( selfID, docString ) &
        bind(C,name='Teuchos_ParameterEntry_setDocString')
    import :: FT_Teuchos_ParameterEntry_ID_t ,c_char
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)              ,intent(in)         ,dimension(*) :: docString
  end subroutine


  !> <BR> Original C++ prototype:
  !! ParameterList& setList( bool isDefault = false, const std::string &docString = "" );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterEntry_setList ( CT_Teuchos_ParameterEntry_ID_t selfID, 
  !!     boolean isDefault, const char docString[] );

  function Teuchos_ParameterEntry_setList ( selfID, isDefault, docString ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_setList')
    import :: FT_Teuchos_ParameterList_ID_t ,FT_Teuchos_ParameterEntry_ID_t ,FT_boolean_t , &
          c_char
    
    type(FT_Teuchos_ParameterList_ID_t)                                   :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)               ,intent(in)   ,value              :: isDefault
    character(kind=c_char)              ,intent(in)         ,dimension(*) :: docString
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> inline T& getValue(T *ptr) const;
  !> <BR> <BR> CTrilinos prototype:
  !! T Teuchos_ParameterEntry_getValue_T ( CT_Teuchos_ParameterEntry_ID_t selfID, T * ptr );

  function Teuchos_ParameterEntry_getValue_double ( selfID, ptr ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_getValue_double')
    import :: c_double ,FT_Teuchos_ParameterEntry_ID_t
    
    real(c_double)                                                        :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                          ,dimension(*) :: ptr
  end function


  !> <BR> Original C++ prototype:
  !! template<typename T> inline T& getValue(T *ptr) const;
  !> <BR> <BR> CTrilinos prototype:
  !! T Teuchos_ParameterEntry_getValue_T ( CT_Teuchos_ParameterEntry_ID_t selfID, T * ptr );

  function Teuchos_ParameterEntry_getValue_int ( selfID, ptr ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_getValue_int')
    import :: c_int ,FT_Teuchos_ParameterEntry_ID_t
    
    integer(c_int)                                                        :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                          ,dimension(*) :: ptr
  end function


  !> <BR> Original C++ prototype:
  !! inline any& getAny(bool activeQry = true);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny ( CT_Teuchos_ParameterEntry_ID_t selfID, 
  !!     boolean activeQry );

  function Teuchos_ParameterEntry_getAny ( selfID, activeQry ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_getAny')
    import :: FT_Teuchos_any_ID_t ,FT_Teuchos_ParameterEntry_ID_t ,FT_boolean_t
    
    type(FT_Teuchos_any_ID_t)                                             :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)               ,intent(in)   ,value              :: activeQry
  end function


  !> <BR> Original C++ prototype:
  !! inline const any& getAny(bool activeQry = true) const;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny_const ( CT_Teuchos_ParameterEntry_ID_t selfID, 
  !!     boolean activeQry );

  function Teuchos_ParameterEntry_getAny_const ( selfID, activeQry ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_getAny_const')
    import :: FT_Teuchos_any_ID_t ,FT_Teuchos_ParameterEntry_ID_t ,FT_boolean_t
    
    type(FT_Teuchos_any_ID_t)                                             :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    integer(FT_boolean_t)               ,intent(in)   ,value              :: activeQry
  end function


  !> <BR> Original C++ prototype:
  !! inline bool isUsed() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterEntry_isUsed ( CT_Teuchos_ParameterEntry_ID_t selfID );

  function Teuchos_ParameterEntry_isUsed ( selfID ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_isUsed')
    import :: FT_boolean_t ,FT_Teuchos_ParameterEntry_ID_t
    
    integer(FT_boolean_t)                                                 :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool isList() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterEntry_isList ( CT_Teuchos_ParameterEntry_ID_t selfID );

  function Teuchos_ParameterEntry_isList ( selfID ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_isList')
    import :: FT_boolean_t ,FT_Teuchos_ParameterEntry_ID_t
    
    integer(FT_boolean_t)                                                 :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! template <typename T> inline bool isType() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterEntry_isType ( CT_Teuchos_ParameterEntry_ID_t selfID );

  function Teuchos_ParameterEntry_isType_double ( selfID ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_isType_double')
    import :: FT_boolean_t ,FT_Teuchos_ParameterEntry_ID_t
    
    integer(FT_boolean_t)                                                 :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! template <typename T> inline bool isType() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterEntry_isType ( CT_Teuchos_ParameterEntry_ID_t selfID );

  function Teuchos_ParameterEntry_isType_int ( selfID ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_isType_int')
    import :: FT_boolean_t ,FT_Teuchos_ParameterEntry_ID_t
    
    integer(FT_boolean_t)                                                 :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! inline bool isDefault() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_ParameterEntry_isDefault ( CT_Teuchos_ParameterEntry_ID_t selfID );

  function Teuchos_ParameterEntry_isDefault ( selfID ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_isDefault')
    import :: FT_boolean_t ,FT_Teuchos_ParameterEntry_ID_t
    
    integer(FT_boolean_t)                                                 :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! inline std::string docString() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Teuchos_ParameterEntry_docString ( CT_Teuchos_ParameterEntry_ID_t selfID );

  function Teuchos_ParameterEntry_docString ( selfID ) result(that) &
        bind(C,name='Teuchos_ParameterEntry_docString')
    import :: c_char ,FT_Teuchos_ParameterEntry_ID_t
    
    character(kind=c_char)                                                :: that
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


!> @name any interface
!! @{

  ! _________________ any interface bodies _________________


  !> <BR> Original C++ prototype:
  !! any();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_any_ID_t Teuchos_any_Create (  );

  function Teuchos_any_Create (  ) result(that) bind(C,name='Teuchos_any_Create')
    import :: FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)                                     :: that
  end function


  !> <BR> Original C++ prototype:
  !! template<typename ValueType> explicit any(const ValueType & value);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_any_ID_t Teuchos_any_Create_ValueType ( ValueType value );

  function Teuchos_any_Create_double ( value ) result(that) &
        bind(C,name='Teuchos_any_Create_double')
    import :: FT_Teuchos_any_ID_t ,c_double
    
    type(FT_Teuchos_any_ID_t)                                     :: that
    real(c_double)              ,intent(in)   ,value              :: value
  end function


  !> <BR> Original C++ prototype:
  !! template<typename ValueType> explicit any(const ValueType & value);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_any_ID_t Teuchos_any_Create_ValueType ( ValueType value );

  function Teuchos_any_Create_int ( value ) result(that) &
        bind(C,name='Teuchos_any_Create_int')
    import :: FT_Teuchos_any_ID_t ,c_int
    
    type(FT_Teuchos_any_ID_t)                                     :: that
    integer(c_int)              ,intent(in)   ,value              :: value
  end function


  !> <BR> Original C++ prototype:
  !! any(const any & other);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_any_ID_t Teuchos_any_Duplicate ( CT_Teuchos_any_ID_t otherID );

  function Teuchos_any_Duplicate ( otherID ) result(that) &
        bind(C,name='Teuchos_any_Duplicate')
    import :: FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)                                     :: that
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: otherID
  end function


  !> <BR> Original C++ prototype:
  !! ~any();
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_any_Destroy ( CT_Teuchos_any_ID_t * selfID );

  subroutine Teuchos_any_Destroy ( selfID ) bind(C,name='Teuchos_any_Destroy')
    import :: FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)                                     :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! any & swap(any & rhs);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Teuchos_any_ID_t Teuchos_any_swap ( CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

  function Teuchos_any_swap ( selfID, rhsID ) result(that) bind(C,name='Teuchos_any_swap')
    import :: FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)                                     :: that
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: rhsID
  end function


  !> <BR> Original C++ prototype:
  !! any & operator=(const any & rhs);
  !> <BR> <BR> CTrilinos prototype:
  !! void Teuchos_any_Assign ( CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

  subroutine Teuchos_any_Assign ( selfID, rhsID ) bind(C,name='Teuchos_any_Assign')
    import :: FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: rhsID
  end subroutine


  !> <BR> Original C++ prototype:
  !! bool empty() const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_any_empty ( CT_Teuchos_any_ID_t selfID );

  function Teuchos_any_empty ( selfID ) result(that) bind(C,name='Teuchos_any_empty')
    import :: FT_boolean_t ,FT_Teuchos_any_ID_t
    
    integer(FT_boolean_t)                                         :: that
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! std::string typeName() const;
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Teuchos_any_typeName ( CT_Teuchos_any_ID_t selfID );

  function Teuchos_any_typeName ( selfID ) result(that) bind(C,name='Teuchos_any_typeName')
    import :: c_char ,FT_Teuchos_any_ID_t
    
    character(kind=c_char)                                        :: that
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! bool same( const any &other ) const;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Teuchos_any_same ( CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t otherID );

  function Teuchos_any_same ( selfID, otherID ) result(that) &
        bind(C,name='Teuchos_any_same')
    import :: FT_boolean_t ,FT_Teuchos_any_ID_t
    
    integer(FT_boolean_t)                                         :: that
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: otherID
  end function


!> @}


  end interface
end module forteuchos

