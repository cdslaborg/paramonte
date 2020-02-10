! Coarray 1D Heat Equation Solver Test: local_field_module
!
! Copyright (c) 2012-2014, Sourcery, Inc.
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
!
module local_field_module
  implicit none
  private
  public :: local_field

  type local_field
    private
    real, allocatable :: values(:)
  contains
    procedure, private :: multiply
    generic :: operator(*)=>multiply
    procedure :: state
    procedure, private :: assign_array
    generic :: assignment(=)=>assign_array
  end type

contains

  pure function multiply(lhs,rhs) result(product_)
    class(local_field), intent(in) :: lhs
    type(local_field) :: product_
    real, intent(in) :: rhs
    product_%values = lhs%values*rhs
  end function

  pure function state(this) result(this_values)
    class(local_field), intent(in) :: this
    real :: this_values(size(this%values))
    this_values = this%values
  end function

  pure subroutine assign_array(lhs,rhs)
    class(local_field), intent(inout) :: lhs
    real, intent(in) :: rhs(:)
    lhs%values = rhs
  end subroutine

end module
