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

module periodic_2nd_order_module
  !note co-object and field modules no longer necessary so they are not used
  use kind_parameters ,only : rkind, ikind
  use ForTrilinos_assertion_utility, only : assert,error_message
  use object_interface, only : object
  use input_file, only:grid_resolution !added input file so that multiple modules/main can access grid_resolution,CZC
  use mpi_module, only :mpi_class !CZC
  use mpi_share, only:mpi_object !CZC
  use shared !CZC
  implicit none
  private
  public :: periodic_2nd_order, initial_field

  type, extends(object) :: periodic_2nd_order
  !   private
     !Make arrays larger than necessary as they must have an explicit intialization
     real(rkind), allocatable :: global_f(:) !object MPI variable for storing function, CZC

  contains
    procedure :: construct
    procedure :: assign_field
    procedure :: add           => add_field
    procedure :: multiply      => multiply_field
    procedure :: multiply_real  !CZC
    procedure :: subtract       !CZC
    procedure :: x             => df_dx
    procedure :: xx            => d2f_dx2
    procedure :: runge_kutta_2nd_step => rk2_dt
    procedure, nopass :: this_image_contains
    procedure :: has_a_zero_at
    procedure :: local_state
    procedure, nopass :: set_time
    procedure, nopass :: get_time
    generic   :: assignment(=) => assign_field
    generic   :: operator(+)   => add
    generic   :: operator(*)   => multiply
    generic   :: operator(*)    => multiply_real!added functionality, taken from field class, CZC
    generic   :: operator(-)    => subtract !added functionality, taken from field class, CZC
    procedure :: output
   !generic :: write=>output ! Fortran 2003 derived-type output
  end type

  real(rkind) ,parameter :: pi=acos(-1._rkind)
  real(rkind), allocatable :: local_grid(:)
  real(rkind) :: time=0.

  abstract interface
    real(rkind) pure function initial_field(x)
      import :: rkind
      real(rkind) ,intent(in) :: x
    end function
  end interface

contains

  pure function local_state(this) result(local_state_vector)
     class(periodic_2nd_order), intent(in) :: this
     real(rkind), allocatable :: local_state_vector(:)
     integer(ikind) :: i
     local_state_vector = this%global_f
  end function

  subroutine output(this,unit,iotype,v_list,iostat,iomsg)
    class(periodic_2nd_order), intent(in) :: this
    integer, intent(in) :: unit ! Unit on which output happens (negative for internal file)
    character(*), intent(in) :: iotype ! Allowable values: ’LISTDIRECTED’,’NAMELIST’, or ’DT’
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    integer(ikind) i
    ! Requires
    call assert(this%user_defined(),error_message("periodic_2nd_order%output recieved unitialized object."))

   do i = 1, local_grid_resolution
      write (unit=unit,iostat=iostat,fmt="(i8,3(f12.4,2x))") &
      (my_id)*local_grid_resolution + i, local_grid(i),time,this%global_f(i) !modified earlier code to work with MPI, CZC
    end do
  end subroutine

  subroutine set_time(time_stamp)
    real(rkind), intent(in) :: time_stamp
    time = time_stamp
  end subroutine

  pure function get_time() result(t)
    real(rkind) :: t
    t = time
  end function

 pure function has_a_zero_at(this, expected_location) result(zero_at_expected_location)
    class(periodic_2nd_order) ,intent(in) :: this
    real(rkind) ,intent(in) :: expected_location
    real(rkind), parameter :: tolerance = 1.0E-06_rkind
    integer :: nearest_grid_point
    logical :: zero_at_expected_location
    ! Requires
    if (this%user_defined()) then
      nearest_grid_point = minloc(abs(local_grid-expected_location),dim=1)
      zero_at_expected_location = merge(.true.,.false., abs(this%global_f(nearest_grid_point)) < tolerance  )
    end if
  end function

  pure function this_image_contains(location) result(within_bounds)
     implicit none
    real(rkind), intent(in) :: location
    logical within_bounds
    within_bounds = merge(.true.,.false., (location>=minval(local_grid) .and. location<=maxval(local_grid)) )
  end function

  subroutine construct (this, initial,num_grid_pts)
     implicit none
    class(periodic_2nd_order), intent(inout) :: this
    procedure(initial_field) ,pointer, intent(in) :: initial
    integer(ikind) ,intent(in) :: num_grid_pts
    integer :: i
    !local variables for storing nearby nodes
    DOUBLE PRECISION left,right !CZC
    DOUBLE PRECISION left_sub,right_sub!CZC
    ! Requires
    call assert(mod(num_grid_pts, num_procs)==0,error_message("periodic_2nd_order%construct: invalid number of grid points."))
    local_grid = grid()
    allocate(this%global_f(local_grid_resolution))
    do concurrent (i=1:local_grid_resolution)
      this%global_f(i) = initial(local_grid(i))
    end do
  ! Ensures

    call mpi_object%oned_message(this%global_f(1:local_grid_resolution),local_grid_resolution,left_sub,right_sub)!communicate with neighbors,

    call this%mark_as_defined
  contains
    pure function grid()
     implicit none
      integer(ikind) :: i
      real(rkind) ,dimension(:) ,allocatable :: grid
      allocate(grid(local_grid_resolution))
      do concurrent (i=1:local_grid_resolution)
        grid(i)  = 2.*pi*(local_grid_resolution*(my_id)+i-1) &
                   /real(num_grid_pts,rkind)
      end do
    end function
  end subroutine

  real(rkind) function rk2_dt(this,nu,num_grid_pts)
     implicit none
    class(periodic_2nd_order) ,intent(in) :: this
    real(rkind) ,intent(in) :: nu
    integer(ikind) ,intent(in) :: num_grid_pts
    real(rkind)             :: dx, CFL, k_max
    ! Requires
    if (this%user_defined()) then
      dx=2.0*pi/num_grid_pts
      k_max=num_grid_pts/2.0_rkind
      CFL=1.0/(1.0-cos(k_max*dx))
      rk2_dt = CFL*dx**2/nu
    end if
  end function

  ! this is the assignment
  subroutine assign_field(lhs,rhs)
    implicit none
    class(periodic_2nd_order) ,intent(inout) :: lhs
    type(periodic_2nd_order) ,intent(in) :: rhs
    DOUBLE PRECISION left,right !CZC
    DOUBLE PRECISION left_sub,right_sub !CZC

    ! Requires
    call assert(rhs%user_defined(),error_message("periodic_2nd_order%copy received undefind RHS."))
    ! update global field
    lhs%global_f = rhs%global_f !CZC
    ! Ensures
    call lhs%mark_as_defined

  end subroutine

  function add_field (this, rhs)
     implicit none
    class(periodic_2nd_order), intent(in) :: this
     class(periodic_2nd_order), intent(in) :: rhs
    type(periodic_2nd_order) :: add_field
   ! Requires
    allocate(add_field%global_f(local_grid_resolution))
    if (rhs%user_defined() .and. this%user_defined()) then
      add_field%global_f(1:local_grid_resolution) = rhs%global_f(1:local_grid_resolution)+this%global_f(1:local_grid_resolution)
      ! Ensures
      call add_field%mark_as_defined
    end if
  end function

  function multiply_field (this, rhs)
     implicit none
    class(periodic_2nd_order), intent(in) :: this, rhs
    type(periodic_2nd_order) :: multiply_field

     ! Requires
    allocate(multiply_field%global_f(local_grid_resolution))
    if (this%user_defined() .and. rhs%user_defined()) then
       multiply_field%global_f(1:local_grid_resolution)=this%global_f(1:local_grid_resolution)*rhs%global_f(1:local_grid_resolution)
      ! Ensures
      call multiply_field%mark_as_defined
    end if
  end function

!New procedure, functionality taken from field, CZC
  function multiply_real(lhs,rhs) result(product_)
    class(periodic_2nd_order) ,intent(in) :: lhs
    real(rkind) ,intent(in)  :: rhs
    type(periodic_2nd_order) :: product_
    ! Requires
    allocate(product_%global_f(local_grid_resolution))
    if (lhs%user_defined()) then
      product_%global_f(1:local_grid_resolution) = lhs%global_f(1:local_grid_resolution) * rhs !multiply array with scalar
      ! Ensures
      call product_%mark_as_defined
    end if
  end function
   !new procedure, functionality taken from field, CZC
  pure function subtract(lhs,rhs) result(difference)
    class(periodic_2nd_order) ,intent(in) :: lhs
    class(periodic_2nd_order) ,intent(in)  :: rhs
    type(periodic_2nd_order) :: difference
    ! Requires
    allocate(difference%global_f(local_grid_resolution))
    if (lhs%user_defined() .and. rhs%user_defined()) then
      difference%global_f(1:local_grid_resolution) = lhs%global_f(1:local_grid_resolution) - rhs%global_f(1:local_grid_resolution) !subtract arrays
      ! Ensures
      call difference%mark_as_defined
    end if
  end function

  function df_dx(this)
     implicit none
    class(periodic_2nd_order), intent(in) :: this
    type(periodic_2nd_order)  :: df_dx
    integer(ikind) :: i,nx
    real(rkind) :: dx, left_image, right_image
    real(rkind), dimension(:), allocatable, save :: tmp_field_array
    ! Requires
    if (this%user_defined()) then

        nx = local_grid_resolution
        if (.not.allocated(tmp_field_array)) allocate(tmp_field_array(nx))
        dx=2.*pi/(real(nx,rkind)*num_procs)

        if (num_procs > 1) then
            call MPI_SENDRECV(this%global_f(1),1,MPI_DOUBLE_PRECISION,left_id,0, &
            right_image,1,MPI_DOUBLE_PRECISION,right_id,0,MPI_COMM_CART,status,ierr)
            call MPI_SENDRECV(this%global_f(local_grid_resolution),1,MPI_DOUBLE_PRECISION,right_id,0, &
            left_image,1,MPI_DOUBLE_PRECISION,left_id,0,MPI_COMM_CART,status,ierr)
        else
            left_image = this%global_f(nx)
            right_image = this%global_f(1)
        end if

      tmp_field_array(1) = &
         0.5*(this%global_f(2)-left_image)/dx

      tmp_field_array(nx) = &
         0.5*(right_image-this%global_f(nx-1))/dx

      do concurrent(i=2:nx-1)
        tmp_field_array(i)=&
          0.5*(this%global_f(i+1)-this%global_f(i-1))/dx
      end do

      df_dx%global_f = tmp_field_array
      ! Ensures
      call df_dx%mark_as_defined
    end if
  end function

  function d2f_dx2(this)
    implicit none
    class(periodic_2nd_order), intent(in) :: this
    type(periodic_2nd_order)  :: d2f_dx2
    integer(ikind) :: i,nx
    real(rkind) :: dx, left_image, right_image
    real(rkind), dimension(:), allocatable, save :: tmp_field_array

    ! Requires
    if (this%user_defined()) then

        nx = local_grid_resolution
        if (.not.allocated(tmp_field_array)) allocate(tmp_field_array(nx))
        dx=2.*pi/(real(nx,rkind)*num_procs)

        if (num_procs > 1) then
            call MPI_SENDRECV(this%global_f(1),1,MPI_DOUBLE_PRECISION,left_id,0, &
            right_image,1,MPI_DOUBLE_PRECISION,right_id,0,MPI_COMM_CART,status,ierr)
            call MPI_SENDRECV(this%global_f(local_grid_resolution),1,MPI_DOUBLE_PRECISION,right_id,0, &
            left_image,1,MPI_DOUBLE_PRECISION,left_id,0,MPI_COMM_CART,status,ierr)
        else
            left_image = this%global_f(nx)
            right_image = this%global_f(1)
        end if

        tmp_field_array(1) = &
         (this%global_f(2)-2.0*this%global_f(1)+left_image)&
         /dx**2

      tmp_field_array(nx) =&
         (right_image-2.0*this%global_f(nx)+this%global_f(nx-1))&
         /dx**2

      do concurrent (i=2:nx-1)
        tmp_field_array(i)=&
          (this%global_f(i+1)-2.0*this%global_f(i)+this%global_f(i-1))&
          /dx**2
      end do

      d2f_dx2%global_f = tmp_field_array
      ! Ensures
      call d2f_dx2%mark_as_defined
    end if
  end function
end module
