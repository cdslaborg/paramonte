! Coarray 1D Heat Equation Solver Test: global_field_module
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
module global_field_module
  use local_field_module, only : local_field
  implicit none
  private
  public :: global_field

  type global_field
    private
    real, allocatable :: values(:)[:]
  contains
    procedure :: set
    procedure :: only_allocate
    generic :: global_field_=>set,only_allocate
    procedure, private :: laplacian
    generic :: operator(.laplacian.) => laplacian
    procedure, private :: add_local_field
    generic :: operator(+) => add_local_field
    procedure, private :: assign_local_field
    generic :: assignment(=) => assign_local_field
    procedure :: state
  end type

  real :: dx
  integer, allocatable :: num_local_points
  integer, parameter:: num_end_points=2
  real :: boundary_vals(num_end_points)

contains

  subroutine only_allocate(this)
    class(global_field), intent(inout) :: this
    if (.not.allocated(num_local_points)) error stop "global_field: no value established for memory allocation yet."
    allocate(this%values(num_local_points)[*]) ! Implicit synchronization point
  end subroutine

  subroutine set(this,internal_values,boundary_values,domain,num_global_points)
    class(global_field), intent(inout) :: this
    integer, intent(in) :: num_global_points
    real, intent(in) :: internal_values,domain(num_end_points),boundary_values(num_end_points)
    if (mod(num_global_points,num_images())/=0) error stop "set: num_global_points not evenly divisible by num_images()"
    if (this_image()==1 .or. this_image()==num_images()) boundary_vals = boundary_values
    if (.not.allocated(num_local_points)) num_local_points=num_global_points/num_images()
    dx=(domain(2)-domain(1))/num_global_points
    allocate(this%values(num_local_points)[*])
    associate(west=>1,east=>2)
      this%values(1) = merge(boundary_values(west),internal_values,this_image()==1)
      this%values(2:num_local_points-1) = internal_values
      this%values(num_local_points) = merge(boundary_values(east),internal_values,this_image()==num_images())
    end associate
    call synchronize()
  end subroutine

  subroutine synchronize()
    if (num_images()>1) then
      associate(me=>this_image())
        if (me==1) then
          sync images(me+1)
        else if (me==num_images()) then
          sync images(me-1)
        else
          sync images([me-1,me+1])
        end if
      end associate
    end if
  end subroutine

  pure function laplacian(rhs) result(laplacian_rhs)
    class(global_field), intent(in) :: rhs
    type(local_field) :: laplacian_rhs
    real :: local_laplacian(num_local_points)
    integer :: i
    associate(N=>num_local_points,me=>this_image())
      if (me==1) then
        local_laplacian(1) = 0.
      else
        local_laplacian(1)=(rhs%values(2)-2.*rhs%values(1)+rhs%values(N)[me-1])/dx**2
      end if
      do concurrent(i=2:N-1)
        local_laplacian(i)=(rhs%values(i+1)-2.*rhs%values(i)+rhs%values(i-1))/dx**2
      end do
      if (me==num_images()) then
        local_laplacian(N) = 0.
      else
        local_laplacian(N)=(rhs%values(1)[me+1]-2.*rhs%values(N)+rhs%values(N-1))/dx**2
      end if
    end associate
    laplacian_rhs = local_laplacian
  end function

  pure function add_local_field(lhs,rhs) result(total)
    class(global_field), intent(in) :: lhs
    type(local_field), intent(in) :: rhs
    type(local_field) :: total
    total = lhs%values + rhs%state()
  end function

  subroutine assign_local_field(lhs,rhs)
    class(global_field), intent(inout) :: lhs
    class(local_field), intent(in) :: rhs
    lhs%values(:) = rhs%state()
    call synchronize()
  end subroutine

  pure function state(this) result(this_values)
    class(global_field), intent(in) :: this
    real :: this_values(size(this%values(:)))
    this_values = this%values
  end function

end module
