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
! Questions? Contact Karla Morris  (knmorri@sandia.gov) or
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************


module FEpetra_BlockMap
#include "ForTrilinos_config.h"
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal ,only : universal
  use ForTrilinos_error ,only : error
  use FEpetra_Comm  ,only : Epetra_Comm
  use iso_c_binding ,only : c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: Epetra_BlockMap ! Expose type/constructors/methods

  type ,extends(universal)      :: Epetra_BlockMap 
    private
    type(FT_Epetra_BlockMap_ID_t) :: BlockMap_id 
  contains
     !Constructors
     procedure ,private :: create_
     procedure ,private :: from_struct_
     generic :: Epetra_BlockMap_ => from_struct_,create_
     !Local/Global ID accessor methods
     !Size and dimension acccessor functions
     procedure         :: NumGlobalElements
     procedure         :: NumMyElements
     procedure         :: IndexBase
     procedure         :: SameAs
     procedure         :: PointSameAs
     procedure         :: MyGlobalElements
     procedure         :: ElementSize_Const
     procedure         :: ElementSize_LID
     generic :: ElementSize=>ElementSize_Const,ElementSize_LID
     !Miscellaneous boolean tests
     procedure         :: LinearMap
     procedure         :: DistributedGlobal
     !Array accessor functions
     !Miscellaneous
     procedure         :: Comm
     !Developers only  -- to be called by developers from other ForTrilinos modules, not by end applications:
     procedure         :: invalidate_id => invalidate_EpetraBlockMap_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraBlockMap
     procedure         :: get_EpetraBlockMap_ID 
     procedure ,nopass :: alias_EpetraBlockMap_ID
     procedure         :: generalize 
  end type

   interface Epetra_BlockMap ! constructors
     !User interface -- constructors for use by end applications:
     module procedure create,duplicate,create_linear,create_arbitrary,create_variable
     !Developers only -- to be called by developers from other ForTrilinos modules, not by end applications:
     module procedure from_struct
   end interface
 
contains
  ! Common type-bound procedures begin here: all ForTrilinos parent classes of the API have procedures analogous to these.
  ! The function from_struct should be called by developers only from other ForTrilinos modules, never by end applications:

  subroutine from_struct_(this,id)
    class(Epetra_BlockMap) ,intent(out) :: this
    type(FT_Epetra_BlockMap_ID_t) ,intent(in) :: id
    this%BlockMap_id = id
    call this%register_self()
  end subroutine

  function from_struct(id) result(new_Epetra_BlockMap)
    type(Epetra_BlockMap) :: new_Epetra_BlockMap
    type(FT_Epetra_BlockMap_ID_t) ,intent(in) :: id
    call new_Epetra_BlockMap%Epetra_BlockMap_(id)
  end function

  ! All additional constructors should take two steps: (1) obtain a struct ID by invoking a procedural binding and then (2) pass
  ! this ID to from_struct to initialize the constructed object's ID component and register the object for reference counting.

  subroutine create_(this,Num_GlobalElements,Element_Size,IndexBase,comm)
    class(Epetra_BlockMap) ,intent(out) :: this 
    integer(c_int) ,intent(in) :: Num_GlobalElements,Element_Size,IndexBase
    class(Epetra_Comm)         :: comm
    call this%Epetra_BlockMap_(Epetra_BlockMap_Create(Num_GlobalElements,Element_Size,IndexBase,comm%get_EpetraComm_ID()))
  end subroutine

  function create(Num_GlobalElements,Element_Size,IndexBase,comm) result(new_Epetra_BlockMap)
    type(Epetra_BlockMap) :: new_Epetra_BlockMap
    integer(c_int) ,intent(in) :: Num_GlobalElements,Element_Size,IndexBase
    class(Epetra_Comm)         :: comm
    call new_Epetra_BlockMap%Epetra_BlockMap_(Num_GlobalElements,Element_Size,IndexBase,comm)
  end function

  type(Epetra_BlockMap) function create_linear(Num_GlobalElements,Num_MyElements,Element_Size,IndexBase,comm)
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) :: Element_Size
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
    type(FT_Epetra_BlockMap_ID_t) :: create_linear_id
    create_linear_id = &
           Epetra_BlockMap_Create_Linear(Num_GlobalElements,Num_MyElements,Element_Size,IndexBase,comm%get_EpetraComm_ID())
    create_linear = from_struct(create_linear_id)
  end function

  type(Epetra_BlockMap) function create_arbitrary(Num_GlobalElements,&
                                                        My_GlobalElements,Element_Size,IndexBase,comm)
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements  
    integer(c_int) ,intent(in) :: Element_Size
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
    type(FT_Epetra_BlockMap_ID_t) :: create_arbitrary_id
    create_arbitrary_id = Epetra_BlockMap_Create_Arbitrary(Num_GlobalElements,size(My_GlobalElements),&
                 My_GlobalElements,Element_Size,IndexBase,comm%get_EpetraComm_ID())
    create_arbitrary = from_struct(create_arbitrary_id)
  end function

  type(Epetra_BlockMap) function create_variable(Num_GlobalElements,&
                                                       My_GlobalElements,Element_SizeList,IndexBase,comm)
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements  
    integer(c_int) ,intent(in) ,dimension(:) :: Element_SizeList    
    integer(c_int) ,intent(in) :: IndexBase
    class(Epetra_Comm)         :: comm
    type(FT_Epetra_BlockMap_ID_t) :: create_variable_id
    create_variable_id = Epetra_BlockMap_Create_Variable(Num_GlobalElements,size(My_GlobalElements),&
                My_GlobalElements,Element_SizeList,IndexBase,comm%get_EpetraComm_ID())
    create_variable = from_struct(create_variable_id)
  end function

  type(Epetra_BlockMap) function duplicate(this)
    type(Epetra_BlockMap) ,intent(in) :: this 
    type(FT_Epetra_BlockMap_ID_t) :: duplicate_id
    duplicate_id = Epetra_BlockMap_Duplicate(this%BlockMap_id)
    duplicate = from_struct(duplicate_id)
  end function

  !----------------- Struct access ---------------------------------------------

  type(FT_Epetra_BlockMap_ID_t) function get_EpetraBlockMap_ID(this)
    class(Epetra_BlockMap) ,intent(in) :: this 
    get_EpetraBlockMap_ID=this%BlockMap_id
  end function

 !----------------- Type casting ---------------------------------------------
  
  type(FT_Epetra_BlockMap_ID_t) function alias_EpetraBlockMap_ID(generic_id)
    use ForTrilinos_table_man, only : CT_Alias
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Epetra_BlockMap_ID
    use iso_c_binding     ,only: c_loc,c_int
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_BlockMap_ID),stat=status)
    ierr=error(status,'FEpetra_BlockMap:alias_EpetraBlockMap_ID')
    call ierr%check_success()
    alias_EpetraBlockMap_ID=degeneralize_EpetraBlockMap(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(Epetra_BlockMap) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%BlockMap_ID))
  end function

 type(FT_Epetra_BlockMap_ID_t) function degeneralize_EpetraBlockMap(generic_id) bind(C)
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_BlockMap_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                   ,value   :: generic_id
    type(FT_Epetra_BlockMap_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraBlockMap = local_ptr
  end function
  
  !_______________Size and dimension acccessor functions____________________
 
  integer(c_int) function NumGlobalElements(this)
    class(Epetra_BlockMap) ,intent(in) :: this
    NumGlobalElements=Epetra_BlockMap_NumGlobalElements(this%BlockMap_id)
  end function 

  integer(c_int) function NumMyElements(this)
    class(Epetra_BlockMap) ,intent(in) :: this
    NumMyElements=Epetra_BlockMap_NumMyElements(this%BlockMap_id)
  end function 

  integer(c_int) function IndexBase(this)
    class(Epetra_BlockMap) ,intent(in) :: this
    IndexBase=Epetra_BlockMap_IndexBase(this%BlockMap_id)
  end function 

  logical function  SameAs(lhs,rhs)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap)        ,intent(in) :: lhs
    class(Epetra_BlockMap)        ,intent(in) :: rhs
    integer(FT_boolean_t) :: SameAs_out
    SameAs_out=Epetra_BlockMap_SameAs(lhs%BlockMap_id,rhs%BlockMap_id)
    if (SameAs_out==FT_FALSE) SameAs=.false.
    if (SameAs_out==FT_TRUE)  SameAs=.true.
  end function SameAs

  logical function  PointSameAs(lhs,rhs)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap)        ,intent(in) :: lhs
    class(Epetra_BlockMap)        ,intent(in) :: rhs
    integer(FT_boolean_t) :: PointSameAs_out
    PointSameAs_out=Epetra_BlockMap_PointSameAs(lhs%BlockMap_id,rhs%BlockMap_id)
    if (PointSameAs_out==FT_FALSE) PointSameAs=.false.
    if (PointSameAs_out==FT_TRUE)  PointSameAs=.true.
  end function PointSameAs

  function MyGlobalElements(this) result(MyGlobalElementsList)
    class(Epetra_BlockMap)     ,intent(in)    :: this
    integer(c_int),dimension(:),allocatable   :: MyGlobalElementsList
    integer(c_int)                            :: junk
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.allocated(MyGlobalElementsList)) then
     allocate(MyGlobalElementsList(this%NumMyElements()),stat=status)
     ierr=error(status,'FEpetra_BlockMap:MyGlobalElements')
     call ierr%check_success()
    endif
    junk=Epetra_BlockMap_MyGlobalElements_Fill(this%BlockMap_id,MyGlobalElementsList)
  end function 

  integer(c_int) function ElementSize_Const(this)
    class(Epetra_BlockMap) ,intent(in) :: this
    ElementSize_Const=Epetra_BlockMap_ElementSize_Const(this%BlockMap_id)
  end function 

  integer(c_int) function ElementSize_LID(this,L_ID)
    class(Epetra_BlockMap) ,intent(in) :: this
    integer(c_int)         ,intent(in) :: L_ID
    integer(c_int)          :: L_ID_c
    L_ID_c=L_ID-FT_Index_OffSet ! To account for Fortran index base 1
    ElementSize_LID=Epetra_BlockMap_ElementSize(this%BlockMap_id,L_ID_c)
  end function 
 
  !_________________Miscellaneous boolean tests____________________________

  logical function LinearMap(this)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap),intent(in) :: this
    integer(FT_boolean_t) :: LinearMap_out
    LinearMap_out=Epetra_BlockMap_LinearMap(this%BlockMap_id)
    if (LinearMap_out==FT_FALSE) LinearMap=.false.
    if (LinearMap_out==FT_TRUE) LinearMap=.true.
  end function

  logical function DistributedGlobal(this)
    use ForTrilinos_enums, only:FT_boolean_t,FT_FALSE,FT_TRUE
    class(Epetra_BlockMap),intent(in) :: this
    integer(FT_boolean_t) :: DistributedGlobal_out
    DistributedGlobal_out=Epetra_BlockMap_DistributedGlobal(this%BlockMap_id)
    if (DistributedGlobal_out==FT_FALSE) DistributedGlobal =.false.
    if (DistributedGlobal_out==FT_TRUE) DistributedGlobal=.true.
  end function

  subroutine Comm(this,communicator)
    use FEpetra_MpiComm
    use FEpetra_SerialComm
    class(Epetra_BlockMap), intent(in) :: this
    class(Epetra_Comm), allocatable :: communicator
#ifdef HAVE_MPI
    type(Epetra_MpiComm), allocatable :: local_Comm
    allocate(Epetra_MpiComm :: local_Comm)
    local_Comm = Epetra_MpiComm(Epetra_BlockMap_Comm(this%BlockMap_id) )
    call move_alloc(local_Comm, communicator)
#else     
    type(Epetra_SerialComm), allocatable :: local_Comm
    allocate(Epetra_SerialComm :: local_Comm)
    local_Comm = Epetra_SerialComm(Epetra_BlockMap_Comm(this%BlockMap_id) )
    call move_alloc(local_Comm,communicator)
#endif
 end subroutine

  !__________ Garbage collection __________________________________________________

  subroutine invalidate_EpetraBlockMap_ID(this)
    class(Epetra_BlockMap),intent(inout) :: this
    this%BlockMap_id%table = FT_Invalid_ID
    this%BlockMap_id%index = FT_Invalid_Index 
    this%BlockMap_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraBlockMap(this)
    class(Epetra_BlockMap),intent(inout) :: this
    call Epetra_BlockMap_Destroy( this%BlockMap_id ) 
  end subroutine

end module 

