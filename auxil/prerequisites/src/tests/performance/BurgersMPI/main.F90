! Copyright (c) 2011, Damian Rouson, Jim Xia, and Xiaofeng Xu.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!     * Neither the names of Damian Rouson, Jim Xia, and Xiaofeng Xu nor the
!       names of any other contributors may be used to endorse or promote products
!       derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL DAMIAN ROUSON, JIM XIA, and XIAOFENG XU BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

module initializer
  use kind_parameters ,only : rkind
  implicit none
contains
  real(rkind) pure function u_initial(x)
    real(rkind) ,intent(in) :: x
    u_initial = 10._rkind*sin(x)
  end function
  real(rkind) pure function zero(x)
    real(rkind) ,intent(in) :: x
    zero = 0.
  end function
end module

program main
  use iso_fortran_env, only : output_unit
  use kind_parameters ,only : rkind
  use periodic_2nd_order_module, only : periodic_2nd_order, initial_field
  use initializer ,only : u_initial,zero
  use input_file ,only : grid_resolution !CZC
  use shared
  use mpi_share, only : mpi_object !CZC
  implicit none
  type(periodic_2nd_order), save :: u,half_uu,u_half
  real(rkind) :: dt,half=0.5,t=0.,t_final=3.08,nu=1.
  integer ,parameter :: base_output_unit=output_unit+10
  integer :: step,iostat, steps, num_steps = 100000
  character(:), allocatable :: iotype ! Allowable values: ’LISTDIRECTED’,’NAMELIST’, or ’DT’
  character(:), allocatable :: iomsg
  integer, allocatable :: v_list(:)
  procedure(initial_field) ,pointer :: initial
  real(rkind), parameter :: time_initial=0.
  real(rkind), allocatable :: u_surface(:,:)
  real(rkind) :: t_1, t_2, t_3

  ! Test parameters
  real(rkind), parameter :: pi=acos(-1._rkind),expected_zero_location=pi

  call mpi_object%mpi_begin !initiate MPI functionality, CZC
  local_grid_resolution=grid_resolution/num_procs !calculate how processors are shared, CZC

#ifdef USING_TAU
  call TAU_PROFILE_SET_NODE(my_id)
#endif
  call cpu_time(t_1)
  initial => u_initial
  call u%construct(initial,grid_resolution)
  initial => zero
  call half_uu%construct(initial,grid_resolution)
  call u_half%construct(initial,grid_resolution)
  call u%set_time(time_initial)
  step = 1
  call u%set_time((step-1)*dt)
  !numerical scheme
  call cpu_time(t_2)
  dt = u%runge_kutta_2nd_step(nu ,grid_resolution)
#ifdef BENCHMARK
  do steps = 1, num_steps
#else
  do while (t<t_final)
#endif
    half_uu = u*u*half
    u_half = u + (u%xx()*nu - half_uu%x())*dt*half ! first substep
    half_uu = u_half*u_half*half
    u  = u + (u_half%xx()*nu - half_uu%x())*dt ! second substep
    t = t + dt
    step = step + 1
  end do
  call cpu_time(t_3)
  if (my_id == 0) print *, t_2 - t_1, t_3 - t_2, t_3 - t_1
  !print *, 'this image = ', my_id, ' f_global = ',u%global_f(1:local_grid_resolution)
  iomsg = "Output result: success."
  !call u%output(70 + my_id,iotype,v_list,iostat,iomsg)
   if (u%this_image_contains(expected_zero_location)) then
     if (.not. u%has_a_zero_at(expected_zero_location)) error stop "Test failed."
     print *,'Test passed.'
   end if
 call mpi_object%mpi_end !end MPI processes

end program
