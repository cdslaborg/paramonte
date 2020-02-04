! Copyright (c) 2012-2016, Sourcery, Inc.
! All rights reserved.
!
! Unit tests for co_min: verify parallel, collective minimum
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!     * Neither the name of Sourcery, Inc., nor the
!       names of any other contributors may be used to endorse or promote products
!       derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL SOURCERY, INC., BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

module co_intrinsics_module
#ifdef USE_EXTENSIONS
  use opencoarrays
#endif
  implicit none

  private
  public :: co_all
  public :: co_product

  interface co_all
    module procedure co_all_logical
  end interface

  interface co_product
    module procedure co_product_c_int,co_product_c_double
  end interface

contains

  subroutine co_all_logical(a)
    logical, intent(inout) :: a(:)
    call co_reduce(a,and)
  contains
    pure function and(lhs,rhs) result(lhs_and_rhs) bind(C,name="and")
      logical, intent(in) :: lhs,rhs
      logical :: lhs_and_rhs
      lhs_and_rhs = lhs .and. rhs
    end function
  end subroutine

  subroutine co_product_c_int(a)
    use iso_c_binding, only : c_int
    integer(c_int), intent(inout) :: a
    call co_reduce(a,product_)
  contains
    pure function product_(lhs,rhs) result(lhs_x_rhs) bind(C,name="product_")
      integer(c_int), intent(in) :: lhs,rhs
      integer(c_int) :: lhs_x_rhs
      lhs_x_rhs = lhs * rhs
    end function
  end subroutine

  subroutine co_product_c_double(a)
    use iso_c_binding, only : c_double
    real(c_double), intent(inout) :: a
    call co_reduce(a,product_)
  contains
    pure function product_(lhs,rhs) result(lhs_x_rhs)
      real(c_double), intent(in) :: lhs,rhs
      real(c_double) :: lhs_x_rhs
      lhs_x_rhs = lhs * rhs
    end function
  end subroutine

end module

program main
  use iso_fortran_env, only : error_unit
  use iso_c_binding, only : c_int,c_double
  use co_intrinsics_module, only : co_all,co_product
#ifdef USE_EXTENSIONS
  use opencoarrays
#endif
  implicit none
  logical :: logical_passes=.false.,c_int_passes=.false.

#ifdef USE_EXTENSIONS
  if (this_image()==1) print *,"Using the extensions from the opencoarrays module."
#endif

  ! Verify that every image has a "true" variable with the value .true.
  verify_co_reduce_logical: block
    logical,dimension(10) :: true=.true.
    sync all
    call co_all(true)
    if (all(true .eqv. .true.)) then
      logical_passes=.true.
    else
      write(error_unit,"(2(a,i2))") "co_reduce fails for logical argument with result (",true,") on image",this_image()
    end if
  end block verify_co_reduce_logical

  ! Verify the product of image number
  verify_co_reduce_c_int: block
    integer(c_int) :: me,i
    me=this_image()
    sync all
    call co_product(me)
    if (me==product(int([(i,i=1,num_images())],c_int))) then
      c_int_passes=.true.
    else
      write(error_unit,"(2(a,i2))") "co_reduce fails integer(c_int) argument with result (",me,") on image",this_image()
    end if
  end block verify_co_reduce_c_int

  ! Verify that this image's tests passed
  if (.not.all([logical_passes,c_int_passes])) error stop

  ! Wait for verification that all images to pass the tests
  sync all
  if (this_image()==1) print *, "Test passed."
end program
