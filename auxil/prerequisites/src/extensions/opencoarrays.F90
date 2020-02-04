! Fortran 2015 feature support for Fortran 2008 compilers
!
! Copyright (c) 2015-2016, Sourcery, Inc.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!     * Neither the name of the Sourcery, Inc., nor the
!       names of its contributors may be used to endorse or promote products
!       derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL SOURCERY, INC., BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
module opencoarrays
#ifdef COMPILER_SUPPORTS_ATOMICS
  use iso_fortran_env, only : atomic_int_kind
#endif
  use iso_c_binding, only : c_int,c_char,c_ptr,c_loc,c_double,c_int32_t,c_ptrdiff_t,c_sizeof,c_bool,c_funloc
  implicit none

#ifndef MPI_WORKING_MODULE
  include 'mpif.h'
#endif

  private
  public :: co_reduce
  public :: co_broadcast
  public :: co_sum
  public :: co_min
  public :: co_max
  public :: this_image
  public :: num_images
  public :: error_stop
  public :: sync_all
  public :: team_number
#ifdef HAVE_MPI
  public :: get_communicator
#endif
#ifdef COMPILER_SUPPORTS_ATOMICS
  public :: event_type
  public :: event_post

  type event_type
    private
    integer(atomic_int_kind), allocatable :: atom[:]
  end type
#endif
#ifdef EXPOSE_INIT_FINALIZE
  public :: caf_init
  public :: caf_finalize
#endif

  ! Generic interface to co_broadcast with implementations for various types, kinds, and ranks
  interface co_reduce
     module procedure co_reduce_c_int,co_reduce_c_double,co_reduce_logical
  end interface

  ! Generic interface to co_broadcast with implementations for various types, kinds, and ranks
  interface co_broadcast
     module procedure co_broadcast_c_int,co_broadcast_c_double,co_broadcast_c_char
  end interface

  ! Generic interface to co_sum with implementations for various types, kinds, and ranks
  interface co_sum
     module procedure co_sum_c_int,co_sum_c_double
  end interface

  ! Generic interface to co_sum with implementations for various types, kinds, and ranks
  interface co_min
     module procedure co_min_c_int,co_min_c_double
  end interface

  ! Generic interface to co_sum with implementations for various types, kinds, and ranks
  interface co_max
     module procedure co_max_c_int,co_max_c_double
  end interface

  abstract interface
     pure function c_int_operator(lhs,rhs) result(lhs_op_rhs)
       import c_int
       integer(c_int), intent(in) :: lhs,rhs
       integer(c_int) :: lhs_op_rhs
     end function
     pure function c_double_operator(lhs,rhs) result(lhs_op_rhs)
       import c_double
       real(c_double), intent(in) :: lhs,rhs
       real(c_double) :: lhs_op_rhs
     end function
     pure function logical_operator(lhs,rhs) result(lhs_op_rhs)
       logical, intent(in) :: lhs,rhs
       logical :: lhs_op_rhs
     end function
  end interface

  interface
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    function team_number(team_type_ptr) result(my_team_number) bind(C,name="_caf_extensions_team_number")
#else
    function team_number(team_type_ptr) result(my_team_number) bind(C,name="_gfortran_caf_team_number")
#endif
       use iso_c_binding, only : c_int,c_ptr
       implicit none
       type(c_ptr), optional :: team_type_ptr
       integer(c_int) :: my_team_number
    end function

#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    function get_communicator(team_type_ptr) result(my_team) bind(C,name="_caf_extensions_get_communicator")
#else
    function get_communicator(team_type_ptr) result(my_team) bind(C,name="_gfortran_caf_get_communicator")
#endif
       use iso_c_binding, only : c_int,c_ptr
       implicit none
       type(c_ptr), optional :: team_type_ptr
       integer(c_int) :: my_team
    end function
  end interface


  ! __________ End Public Interface _____________


  ! __________ Begin Private Implementation _____

  ! Enumeration from ../libcaf-gfortran-descriptor.h:
  !  enum
  !{ BT_UNKNOWN = 0, BT_INTEGER, BT_LOGICAL, BT_REAL, BT_COMPLEX,
  !  BT_DERIVED, BT_CHARACTER, BT_CLASS, BT_PROCEDURE, BT_HOLLERITH, BT_VOID,
  !  BT_ASSUMED
  !};

  enum ,bind(C)
    enumerator :: &
     BT_UNKNOWN = 0, BT_INTEGER, BT_LOGICAL, BT_REAL, BT_COMPLEX, &
     BT_DERIVED, BT_CHARACTER, BT_CLASS, BT_PROCEDURE, BT_HOLLERITH, BT_VOID, &
     BT_ASSUMED
  end enum

  ! Type definition from ../libcaf-gfortran-descriptor.h:
  !typedef struct descriptor_dimension
  !{
  !  ptrdiff_t _stride;
  !  ptrdiff_t lower_bound;
  !  ptrdiff_t _ubound;
  !}
  !descriptor_dimension;

  ! Fortran derived type interoperable with like-named C type:
  type, bind(C) :: descriptor_dimension
     integer(c_ptrdiff_t) :: stride
     integer(c_ptrdiff_t) :: lower_bound
     integer(c_ptrdiff_t) :: ubound_
  end type

  ! Type definition from ../libcaf-gfortran-descriptor.h:
  !typedef struct gfc_descriptor_t {
  !  void *base_addr;
  !  size_t offset;
  !  ptrdiff_t dtype;
  !  descriptor_dimension dim[];
  !} gfc_descriptor_t;

  integer, parameter :: max_dimensions=15

  ! Fortran derived type interoperable with like-named C type:
  type, bind(C) :: gfc_descriptor_t
    type(c_ptr) :: base_addr
    integer(c_ptrdiff_t) :: offset
    integer(c_ptrdiff_t) :: dtype
    type(descriptor_dimension) :: dim_(max_dimensions)
  end type

  ! C comment and source from ../libcaf.h
  ! /* When there is a vector subscript in this dimension, nvec == 0, otherwise,
  ! lower_bound, upper_bound, stride contains the bounds relative to the declared
  ! bounds; kind denotes the integer kind of the elements of vector[].  */
  ! type, bind(C) :: caf_vector_t {
  !   size_t nvec;
  !   union {
  !     struct {
  !       void *vector;
  !       int kind;
  !     } v;
  !     struct {
  !       ptrdiff_t lower_bound, upper_bound, stride;
  !     } triplet;
  !   } u;
  ! }
  ! caf_vector_t;

  type, bind(C) :: v_t
    type(c_ptr) :: vector
    integer(c_int) :: kind_
  end type

  type, bind(C) :: triplet_t
    integer(c_ptrdiff_t) :: lower_bound, upper_bound, stride
  end type

  type, bind(C) :: u_t
     type(v_t) :: v
     type(triplet_t) :: triplet
  end type

  type, bind(C) :: caf_vector_t
    integer(c_ptrdiff_t) :: nvec
    type(u_t) :: u
  end type

  ! --------------------

  integer(c_int), save, volatile, bind(C,name="CAF_COMM_WORLD") :: CAF_COMM_WORLD
  integer(c_int32_t), parameter  :: bytes_per_word=4_c_int32_t

  interface gfc_descriptor
    module procedure gfc_descriptor_c_int,gfc_descriptor_c_double,gfc_descriptor_logical
  end interface

  ! Bindings for OpenCoarrays C procedures
  interface

    ! C function signature from ../mpi/mpi_caf.c:
    ! void
    ! PREFIX (init) (int *argc, char ***argv)
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    subroutine caf_init(argc,argv) bind(C,name="_caf_extensions_init")
#else
    subroutine caf_init(argc,argv) bind(C,name="_gfortran_caf_init")
#endif
      import :: c_int,c_ptr
      type(c_ptr), value :: argc
      type(c_ptr), value :: argv
    end subroutine

    ! C function signature from ../mpi/mpi_caf.c:
    ! void
    ! PREFIX (finalize) (void)
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    subroutine caf_finalize() bind(C,name="_caf_extensions_finalize")
#else
    subroutine caf_finalize() bind(C,name="_gfortran_caf_finalize")
#endif
    end subroutine

    ! C function signature from ../mpi/mpi_caf.c:
    ! void
    ! PREFIX (co_min) (gfc_descriptor_t *a, int result_image, int *stat, char *errmsg,
    !                  int src_len, int errmsg_len)
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    subroutine opencoarrays_co_min(a,result_image, stat, errmsg, unused , errmsg_len) bind(C,name="_caf_extensions_co_min")
#else
    subroutine opencoarrays_co_min(a,result_image, stat, errmsg, unused, errmsg_len) bind(C,name="_gfortran_caf_co_min")
#endif
      import :: c_int,c_char,c_ptr
      type(c_ptr), intent(in), value :: a
      integer(c_int), intent(in),  value :: result_image,errmsg_len,unused
      integer(c_int), intent(out), optional, volatile :: stat
      character(kind=c_char), intent(out), optional, volatile :: errmsg(*)
    end subroutine

    ! C function signature from ../mpi/mpi_caf.c
    ! void
    ! PREFIX (co_reduce) (gfc_descriptor_t *a, void *(*opr) (void *, void *), int opr_flags,
    !                     int result_image, int *stat, char *errmsg, int a_len, int errmsg_len)
    subroutine opencoarrays_co_reduce(a, opr, opr_flags, result_image, stat, errmsg, a_len, errmsg_len) &
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
      bind(C,name="_caf_extensions_co_reduce")
#else
      bind(C,name="_gfortran_caf_co_reduce")
#endif
      use iso_c_binding, only : c_ptr,c_funptr,c_int,c_char
      type(c_ptr), intent(in), value :: a
      type(c_funptr), intent(in), value :: opr
      integer(c_int), intent(in), value :: opr_flags,result_image
      integer(c_int), intent(out) :: stat
      character(kind=c_char), intent(out), optional, volatile :: errmsg(*)
      integer(c_int), intent(in), value :: a_len
      integer(c_int), intent(in), value :: errmsg_len
    end subroutine

    ! C function signature from ../mpi/mpi_caf.c:
    ! void
    ! PREFIX (get) (caf_token_t token, size_t offset,
    !              int image_index,
    !              gfc_descriptor_t *src ,
    !              caf_vector_t *src_vector __attribute__ ((unused)),
    !              gfc_descriptor_t *dest, int src_kind, int dst_kind,
    !              bool mrt)
    subroutine opencoarrays_get(token, offset, image_index_, src, src_vector_unused, dest, src_kind, dst_kind, mrt) &
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
      bind(C,name="_caf_extensions_get")
#else
      bind(C,name="_gfortran_caf_get")
#endif
      import c_ptrdiff_t,c_int,gfc_descriptor_t,c_bool,caf_vector_t,c_ptr
      type(c_ptr), value :: token
      integer(c_ptrdiff_t), value :: offset
      integer(c_int), value :: image_index_
      type(gfc_descriptor_t) :: src
      type(caf_vector_t) :: src_vector_unused
      type(gfc_descriptor_t) :: dest
      integer(c_int), value :: src_kind
      integer(c_int), value :: dst_kind
      logical(c_bool), value :: mrt
    end subroutine

    ! C function signature from ../mpi/mpi_caf.c:
    ! void
    ! PREFIX (co_max) (gfc_descriptor_t *a, int result_image, int *stat,
    !                  char *errmsg, int src_len, int errmsg_len)
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    subroutine opencoarrays_co_max(a,result_image, stat, errmsg, unused, errmsg_len) bind(C,name="_caf_extensions_co_max")
#else
    subroutine opencoarrays_co_max(a,result_image, stat, errmsg, unused, errmsg_len) bind(C,name="_gfortran_caf_co_max")
#endif
      import :: c_int,c_char,c_ptr
      type(c_ptr), intent(in), value :: a
      integer(c_int), intent(in),  value :: result_image,errmsg_len,unused
      integer(c_int), intent(out), optional :: stat
      character(kind=c_char), intent(out), optional :: errmsg(*)
    end subroutine

    ! C function signature from ../mpi/mpi_caf.c:
    ! void
    ! PREFIX (co_sum) (gfc_descriptor_t *a, int result_image, int *stat, char *errmsg,
    !                  int errmsg_len)
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    subroutine opencoarrays_co_sum(a, result_image, stat, errmsg, errmsg_len) bind(C,name="_caf_extensions_co_sum")
#else
    subroutine opencoarrays_co_sum(a, result_image, stat, errmsg, errmsg_len) bind(C,name="_gfortran_caf_co_sum")
#endif
      import :: c_int,c_char,c_ptr
      type(c_ptr), intent(in), value :: a
      integer(c_int), intent(in),  value :: result_image,errmsg_len
      integer(c_int), intent(out), optional :: stat
      character(kind=c_char), intent(out), optional :: errmsg(*)
    end subroutine

    ! C function signature from ../mpi/mpi_caf.c
    ! void
    ! PREFIX (co_broadcast) (gfc_descriptor_t *a, int source_image, int *stat, char *errmsg,
    !                  int errmsg_len)
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    subroutine opencoarrays_co_broadcast(a,source_image, stat, errmsg, errmsg_len) bind(C,name="_caf_extensions_co_broadcast")
#else
    subroutine opencoarrays_co_broadcast(a,source_image, stat, errmsg, errmsg_len) bind(C,name="_gfortran_caf_co_broadcast")
#endif
      import :: c_int,c_char,c_ptr
      type(c_ptr), intent(in), value :: a
      integer(c_int), intent(in),  value :: source_image,errmsg_len
      integer(c_int), intent(out), optional :: stat
      character(kind=c_char), intent(out), optional :: errmsg(*)
    end subroutine

    ! C function signature from ../mpi/mpi_caf.c:
    ! int PREFIX (this_image) (int);
    function opencoarrays_this_image(coarray) bind(C,name="_gfortran_caf_this_image") result(image_num)
      import :: c_int
      integer(c_int), value, intent(in) :: coarray
      integer(c_int)  :: image_num
    end function

    ! C function signature from ../mpi/mpi_caf.c:
    ! int PREFIX (num_images) (int, int);
    function opencoarrays_num_images(coarray,dim_) bind(C,name="_gfortran_caf_num_images") result(num_images_)
      import :: c_int
      integer(c_int), value, intent(in) :: coarray,dim_
      integer(c_int) :: num_images_
    end function

    ! C function signature from ../mpi_caf.c
    ! void PREFIX (error_stop) (int32_t) __attribute__ ((noreturn));
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    subroutine opencoarrays_error_stop(stop_code) bind(C,name="_caf_extensions_error_stop")
#else
    subroutine opencoarrays_error_stop(stop_code) bind(C,name="_gfortran_caf_error_stop")
#endif
      import :: c_int32_t
      integer(c_int32_t), value, intent(in) :: stop_code
    end subroutine

    ! C function signature from ../mpi_caf.c
    ! void PREFIX (sync_all) (int *, char *, int);
#ifdef COMPILER_SUPPORTS_CAF_INTRINSICS
    subroutine opencoarrays_sync_all(stat,errmsg,unused) bind(C,name="_caf_extensions_sync_all")
#else
    subroutine opencoarrays_sync_all(stat,errmsg,unused) bind(C,name="_gfortran_caf_sync_all")
#endif
      import :: c_int,c_char
      integer(c_int), intent(out) :: stat,unused
      character(c_char), intent(out) :: errmsg(*)
    end subroutine

  end interface


contains

  ! __________ Descriptor constructors for for each supported type and kind ____________
  ! ____________________________________________________________________________________

  function my_dtype(type_,kind_,rank_) result(dtype_)
    integer,parameter :: GFC_DTYPE_SIZE_SHIFT = 8, GFC_DTYPE_TYPE_SHIFT=3
    integer(c_int32_t), intent(in) :: type_,kind_,rank_
    integer(c_int32_t) :: dtype_

    ! SIZE Type Rank
    ! 0000  000  000

    ! Rank is represented in the 3 least significant bits
    dtype_ = ior(0_c_int32_t,rank_)
    ! The next three bits represent the type id as expressed in libcaf-gfortran-descriptor.h
    dtype_ = ior(dtype_,ishft(type_,GFC_DTYPE_TYPE_SHIFT))
    ! The most significant bits represent the size of a the type (single or double precision).
    ! We can express the precision in terms of 32-bit words: 1 for single, 2 for double.
    dtype_ = ior(dtype_,ishft(kind_,GFC_DTYPE_SIZE_SHIFT))

  end function

  function gfc_descriptor_c_int(a) result(a_descriptor)
    integer(c_int), intent(in), target, contiguous :: a(..)
    type(gfc_descriptor_t) :: a_descriptor
    integer(c_int), parameter :: unit_stride=1,scalar_offset=-1
    integer(c_int) :: i

    a_descriptor%dtype = my_dtype(type_=BT_INTEGER,kind_=int(c_sizeof(a)/bytes_per_word,c_int32_t),rank_=rank(a))
    a_descriptor%offset = scalar_offset
    a_descriptor%base_addr = c_loc(a) ! data
    do concurrent(i=1:rank(a))
      a_descriptor%dim_(i)%stride  = unit_stride
      a_descriptor%dim_(i)%lower_bound = lbound(a,i)
      a_descriptor%dim_(i)%ubound_ = ubound(a,i)
    end do

  end function

  function gfc_descriptor_logical(a) result(a_descriptor)
    logical, intent(in), target, contiguous :: a(..)
    type(gfc_descriptor_t) :: a_descriptor
    integer(c_int), parameter :: unit_stride=1,scalar_offset=-1,words=1
    integer(c_int) :: i

    a_descriptor%dtype = my_dtype(type_=BT_LOGICAL,kind_=words,rank_=rank(a))
    a_descriptor%offset = scalar_offset
    a_descriptor%base_addr = c_loc(a) ! data
    do concurrent(i=1:rank(a))
      a_descriptor%dim_(i)%stride  = unit_stride
      a_descriptor%dim_(i)%lower_bound = lbound(a,i)
      a_descriptor%dim_(i)%ubound_ = ubound(a,i)
    end do

  end function

  function gfc_descriptor_c_double(a) result(a_descriptor)
    real(c_double), intent(in), target, contiguous :: a(..)
    type(gfc_descriptor_t) :: a_descriptor
    integer(c_int), parameter :: unit_stride=1,scalar_offset=-1
    integer(c_int) :: i

    a_descriptor%dtype = my_dtype(type_=BT_REAL,kind_=int(c_sizeof(a)/bytes_per_word,c_int32_t),rank_=rank(a))
    a_descriptor%offset = scalar_offset
    a_descriptor%base_addr = c_loc(a) ! data
    do concurrent(i=1:rank(a))
      a_descriptor%dim_(i)%stride  = unit_stride
      a_descriptor%dim_(i)%lower_bound = lbound(a,i)
      a_descriptor%dim_(i)%ubound_ = ubound(a,i)
    end do

  end function

 ! This version should work for any rank but causes an ICE with gfortran 4.9.2
 !
 !function gfc_descriptor_c_char(a) result(a_descriptor)
 !  character(c_char), intent(in), target, contiguous :: a(..)
 !  type(gfc_descriptor_t) :: a_descriptor
 !  integer(c_int), parameter :: unit_stride=1,scalar_offset=-1
 !  integer(c_int) :: i

 !  a_descriptor%dtype = my_dtype(type_=BT_CHARACTER,kind_=int(c_sizeof(a)/bytes_per_word,c_int32_t),rank_=rank(a))
 !  a_descriptor%offset = scalar_offset
 !  a_descriptor%base_addr = c_loc(a) ! data
 !  do concurrent(i=1:rank(a))
 !    a_descriptor%dim_(i)%stride  = unit_stride
 !    a_descriptor%dim_(i)%lower_bound = lbound(a,i)
 !    a_descriptor%dim_(i)%ubound_ = ubound(a,i)
 !  end do

 !end function

  ! ______ Assumed-rank co_reduce wrappers for each supported type and kind _________
  ! _________________________________________________________________________________

    subroutine co_reduce_c_int(a, opr, result_image, stat, errmsg)
      ! Dummy variables
      integer(c_int), intent(inout), volatile, contiguous :: a(..)
      procedure(c_int_operator) :: opr
      integer(c_int), intent(in), optional :: result_image
      integer(c_int), intent(out), optional, volatile :: stat
      character(kind=c_char), intent(out), optional, volatile :: errmsg(*)
      ! Local variables
      integer(c_int), volatile :: opr_flags_unused,a_len_unused,errmsg_len
      type(gfc_descriptor_t), target :: a_descriptor
      integer(c_int), parameter :: default_result_image=0
      integer(c_int) :: result_image_

      result_image_ = merge(result_image,default_result_image,present(result_image))
      a_descriptor = gfc_descriptor(a)
      call opencoarrays_co_reduce( &
        c_loc(a_descriptor), c_funloc(opr), opr_flags_unused, result_image_, stat, errmsg, a_len_unused, errmsg_len &
      )
    end subroutine

    subroutine co_reduce_logical(a, opr, result_image, stat, errmsg)
      ! Dummy variables
      logical, intent(inout), volatile, contiguous :: a(..)
      procedure(logical_operator) :: opr
      integer(c_int), intent(in), optional :: result_image
      integer(c_int), intent(out),optional , volatile :: stat
      character(kind=c_char), intent(out), optional, volatile :: errmsg(*)
      ! Local variables
      integer(c_int), volatile :: opr_flags_unused,a_len_unused,errmsg_len
      type(gfc_descriptor_t), target :: a_descriptor
      integer(c_int), parameter :: default_result_image=0
      integer(c_int) :: result_image_

      result_image_ = merge(result_image,default_result_image,present(result_image))
      a_descriptor = gfc_descriptor(a)
      call opencoarrays_co_reduce( &
        c_loc(a_descriptor), c_funloc(opr), opr_flags_unused, result_image_, stat, errmsg, a_len_unused, errmsg_len &
      )
    end subroutine

    subroutine co_reduce_c_double(a, opr, result_image, stat, errmsg)
      ! Dummy variables
      real(c_double), intent(inout), volatile, contiguous :: a(..)
      procedure(c_double_operator) :: opr
      integer(c_int), intent(in), optional :: result_image
      integer(c_int), intent(out), optional, volatile :: stat
      character(kind=c_char), intent(out), optional, volatile :: errmsg(*)
      ! Local variables
      integer(c_int), volatile :: opr_flags_unused,a_len_unused,errmsg_len
      type(gfc_descriptor_t), target :: a_descriptor
      integer(c_int), parameter :: default_result_image=0
      integer(c_int) :: result_image_

      result_image_ = merge(result_image,default_result_image,present(result_image))
      a_descriptor = gfc_descriptor(a)
      call opencoarrays_co_reduce( &
        c_loc(a_descriptor), c_funloc(opr), opr_flags_unused, result_image_, stat, errmsg, a_len_unused, errmsg_len &
      )
    end subroutine

  ! ______ Assumed-rank co_broadcast wrappers for each supported type and kind _________
  ! ____________________________________________________________________________________

  ! This provisional implementation incurs some overhead by converting the character argument
  ! to an integer(c_int) array, invoking co_broadcast_c_int and then convering the received
  ! message back from the integer(c_int) array to a character variable.
  !
  ! Replace this implementation with one that avoids the conversions and the associated copies
  ! once the compiler provides support for co_broadcast with scalar arguments.
  !
  subroutine co_broadcast_c_char(a,source_image,stat,errmsg)
    character(kind=c_char,len=*), intent(inout), volatile, target :: a
    integer(c_int), intent(in), optional :: source_image
    integer(c_int), intent(out), optional:: stat
    character(kind=1,len=*), intent(out), optional :: errmsg
    ! Local variables and constants:
    integer(c_int), allocatable :: a_cast_to_integer_array(:)
    character(len=len(a)) :: dummy

    ! Convert "a" to an integer(c_int) array where each 32-bit integer element holds four 1-byte characters
    a_cast_to_integer_array = transfer(a,[0_c_int])
    ! Broadcast the integer(c_int) array
    call co_broadcast_c_int(a_cast_to_integer_array,source_image, stat, errmsg)
    ! Recover the characters from the broadcasted integer(c_int) array
    dummy = ' '
    a = transfer(a_cast_to_integer_array,dummy)

  end subroutine

  subroutine co_broadcast_c_double(a,source_image,stat,errmsg)
    real(c_double), intent(inout), volatile, target, contiguous :: a(..)
    integer(c_int), intent(in), optional :: source_image
    integer(c_int), intent(out), optional:: stat
    character(kind=1,len=*), intent(out), optional :: errmsg
    ! Local variables and constants
    integer(c_int), parameter :: default_source_image=0
    integer(c_int) :: source_image_ ! Local replacement for the corresponding intent(in) dummy argument
    type(gfc_descriptor_t), target :: a_descriptor

    source_image_ = merge(source_image,default_source_image,present(source_image))
    a_descriptor = gfc_descriptor(a)
    call opencoarrays_co_broadcast(c_loc(a_descriptor),source_image_, stat, errmsg, len(errmsg))

  end subroutine

  subroutine co_broadcast_c_int(a,source_image,stat,errmsg)
    integer(c_int), intent(inout), volatile, target, contiguous :: a(..)
    integer(c_int), intent(in), optional :: source_image
    integer(c_int), intent(out), optional:: stat
    character(kind=1,len=*), intent(out), optional :: errmsg
    ! Local variables and constants:
    integer(c_int), parameter :: default_source_image=0
    integer(c_int) :: source_image_ ! Local replacement for the corresponding intent(in) dummy argument
    type(gfc_descriptor_t), target :: a_descriptor

    source_image_ = merge(source_image,default_source_image,present(source_image))
    a_descriptor = gfc_descriptor(a)
    call opencoarrays_co_broadcast(c_loc(a_descriptor),source_image_, stat, errmsg, len(errmsg))

  end subroutine

  ! ________ Assumed-rank get wrappers for each supported type and kind ________________
  ! ________ (Incomplete, private and unsupported) _____________________________________
  ! ____________________________________________________________________________________

  subroutine get_c_int(src , dest, image_index_, offset, mrt)
    use iso_fortran_env, only : error_unit
    ! Dummy arguments:
    integer(c_int), intent(in), target, contiguous :: src(..)
    integer(c_int), intent(out), target, contiguous, volatile :: dest(..)
    integer(c_int), intent(in) :: image_index_
    integer(c_ptrdiff_t), intent(in) :: offset
    logical(c_bool), intent(in) :: mrt
    ! Local variables:
    type(gfc_descriptor_t), target, volatile ::  dest_descriptor
    type(gfc_descriptor_t), target :: src_descriptor
    type(caf_vector_t) :: src_vector_unused
    type(c_ptr) :: token

    write(error_unit,*) "Remote access of coarrays not yet supported"
    call error_stop

    src_descriptor = gfc_descriptor(src)
    dest_descriptor = gfc_descriptor(dest)
    call opencoarrays_get( &
     token, offset, image_index_, src_descriptor, src_vector_unused, dest_descriptor, kind(src), kind(dest), mrt &
    )
  end subroutine

  ! ________ Assumed-rank co_min wrappers for each supported type and kind _____________
  ! ____________________________________________________________________________________

  subroutine co_min_c_int(a,result_image,stat,errmsg)
    integer(c_int), intent(inout), volatile, target, contiguous :: a(..)
    integer(c_int), intent(in), optional :: result_image
    integer(c_int), intent(out), optional:: stat
    character(kind=1,len=*), intent(out), optional :: errmsg
    ! Local variables and constants:
    type(gfc_descriptor_t), target :: a_descriptor
    integer(c_int), parameter :: default_result_image=0
    integer(c_int) :: unused, result_image_ ! Local replacement for the corresponding intent(in) dummy argument

    a_descriptor = gfc_descriptor(a)
    result_image_ = merge(result_image,default_result_image,present(result_image))
    call opencoarrays_co_min(c_loc(a_descriptor),result_image_, stat, errmsg, unused, len(errmsg))

  end subroutine

  subroutine co_min_c_double(a,result_image,stat,errmsg)
    real(c_double), intent(inout), volatile, target, contiguous :: a(..)
    integer(c_int), intent(in), optional :: result_image
    integer(c_int), intent(out), optional:: stat
    character(kind=1,len=*), intent(out), optional :: errmsg
    ! Local variables and constants:
    type(gfc_descriptor_t), target :: a_descriptor
    integer(c_int), parameter :: default_result_image=0
    integer(c_int) :: unused, result_image_ ! Local replacement for the corresponding intent(in) dummy argument

    a_descriptor = gfc_descriptor(a)
    result_image_ = merge(result_image,default_result_image,present(result_image))
    call opencoarrays_co_min(c_loc(a_descriptor),result_image_, stat, errmsg, unused, len(errmsg))

  end subroutine


  ! ________ Assumed-rank co_max wrappers for each supported type and kind _____________
  ! ____________________________________________________________________________________

  subroutine co_max_c_int(a,result_image,stat,errmsg)
    integer(c_int), intent(inout), volatile, target, contiguous :: a(..)
    integer(c_int), intent(in), optional :: result_image
    integer(c_int), intent(out), optional:: stat
    character(kind=1,len=*), intent(out), optional :: errmsg
    ! Local variables and constants:
    type(gfc_descriptor_t), target :: a_descriptor
    integer(c_int), parameter :: default_result_image=0
    integer(c_int) :: unused, result_image_ ! Local replacement for the corresponding intent(in) dummy argument

    a_descriptor = gfc_descriptor(a)
    result_image_ = merge(result_image,default_result_image,present(result_image))
    call opencoarrays_co_max(c_loc(a_descriptor),result_image_, stat, errmsg, unused, len(errmsg))

  end subroutine

  subroutine co_max_c_double(a,result_image,stat,errmsg)
    real(c_double), intent(inout), volatile, target, contiguous :: a(..)
    integer(c_int), intent(in), optional :: result_image
    integer(c_int), intent(out), optional:: stat
    character(kind=1,len=*), intent(out), optional :: errmsg
    ! Local variables and constants:
    type(gfc_descriptor_t), target :: a_descriptor
    integer(c_int), parameter :: default_result_image=0
    integer(c_int) :: unused, result_image_ ! Local replacement for the corresponding intent(in) dummy argument

    a_descriptor = gfc_descriptor(a)
    result_image_ = merge(result_image,default_result_image,present(result_image))
    call opencoarrays_co_max(c_loc(a_descriptor),result_image_, stat, errmsg, unused, len(errmsg))

  end subroutine

  ! ________ Assumed-rank co_sum wrappers for each supported type and kind _____________
  ! ____________________________________________________________________________________

  subroutine co_sum_c_double(a,result_image,stat,errmsg)
    real(c_double), intent(inout), volatile, target, contiguous :: a(..)
    integer(c_int), intent(in), optional :: result_image
    integer(c_int), intent(out), optional:: stat
    character(kind=1,len=*), intent(out), optional :: errmsg
    ! Local variables and constants:
    type(gfc_descriptor_t), target :: a_descriptor
    integer(c_int), parameter :: default_result_image=0
    integer(c_int) :: result_image_ ! Local replacement for the corresponding intent(in) dummy argument

    a_descriptor = gfc_descriptor(a)
    result_image_ = merge(result_image,default_result_image,present(result_image))
    call opencoarrays_co_sum(c_loc(a_descriptor),result_image_, stat, errmsg, len(errmsg))

  end subroutine

  subroutine co_sum_c_int(a,result_image,stat,errmsg)
    integer(c_int), intent(inout), volatile, target, contiguous :: a(..)
    integer(c_int), intent(in), optional :: result_image
    integer(c_int), intent(out), optional:: stat
    character(kind=1,len=*), intent(out), optional :: errmsg
    ! Local variables and constants:
    integer(c_int), parameter :: default_result_image=0
    type(gfc_descriptor_t), target :: a_descriptor
    integer(c_int) :: result_image_ ! Local replacement for the corresponding intent(in) dummy argument

    a_descriptor = gfc_descriptor(a)
    result_image_ = merge(result_image,default_result_image,present(result_image))
    call opencoarrays_co_sum(c_loc(a_descriptor),result_image_, stat, errmsg, len(errmsg))

  end subroutine

  ! Return the image number (MPI rank + 1)
  function this_image()  result(image_num)
#ifdef MPI_WORKING_MODULE
    use mpi, only : MPI_Comm_rank
#endif
    integer(c_int) :: image_num,ierr
   !image_num = opencoarrays_this_image(unused)
    call MPI_Comm_rank(CAF_COMM_WORLD,image_num,ierr)
    if (ierr/=0) call error_stop
    image_num = image_num + 1
  end function

  ! Return the total number of images
  function num_images()  result(num_images_)
#ifdef MPI_WORKING_MODULE
    use mpi, only : MPI_Comm_size
#endif
    integer(c_int) :: num_images_,ierr
   !num_images_ = opencoarrays_num_images(unused_coarray,unused_scalar)
    call MPI_Comm_size(CAF_COMM_WORLD,num_images_,ierr)
    if (ierr/=0) call error_stop
  end function

  ! Halt the execution of all images
  subroutine error_stop(stop_code)
    integer(c_int32_t), intent(in), optional :: stop_code
    integer(c_int32_t), parameter :: default_code=-1_c_int32_t
    integer(c_int32_t) :: code
    code = merge(stop_code,default_code,present(stop_code))
    call opencoarrays_error_stop(code)
  end subroutine

  ! Impose a global execution barrier
  subroutine sync_all(stat,errmsg,unused)
    integer(c_int), intent(out), optional :: stat,unused
    character(c_char), intent(out), optional :: errmsg
    call opencoarrays_sync_all(stat,errmsg,unused)
  end subroutine

#ifdef COMPILER_SUPPORTS_ATOMICS
   ! Proposed Fortran 2015 event_post procedure
   subroutine event_post(this)
     class(event_type), intent(inout) ::  this
     if (.not.allocated(this%atom)) this%atom=0
     call atomic_define ( this%atom, this%atom + 1_atomic_int_kind )
   end subroutine
#endif

end module
