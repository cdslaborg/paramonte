! Coarray 1D Heat Equation Solver Test: main
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

program main
  use IEEE_arithmetic, only : IEEE_is_NaN
  use global_field_module, only : global_field
  implicit none
  type(global_field) :: T,laplacian_T,T_half
  real, parameter :: alpha=1.,dt=0.0001,final_time=1.,tolerance=1.E-3
  real :: time=0.
  call T%global_field_(internal_values=0.,boundary_values=[1.,0.],domain=[0.,1.],num_global_points=32)
  call T_half%global_field_()
  do while(time<final_time)
    T_half = T + (.laplacian.T)*(alpha*dt/2.)
    T      = T + (.laplacian.T_half)*(alpha*dt)
    time = time + dt
  end do
  call laplacian_T%global_field_()
  laplacian_T = .laplacian.T
  block 
    real, allocatable :: residual(:)
    residual = laplacian_T%state()
    if ( any(residual>tolerance) .or. any(IEEE_is_NaN(residual)) .or. any(residual<0) ) error stop "Test failed."
  end block
  if (this_image()==1) print *,"Test passed."
end program
