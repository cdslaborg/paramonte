module global_field_module
  use iso_fortran_env, only : real64,int64
  use co_object_interface, only : co_object
  use ForTrilinos_assertion_utility, only : assert,error_message
  use local_field_module, only : local_field
  implicit none
  private
  public :: global_field,initial_condition

  type, extends(co_object) :: global_field
    private
    real(real64), allocatable :: values(:)[:]
  contains
    procedure :: set
    procedure :: state
    procedure :: x
    procedure :: xx
    procedure, nopass :: grid_spacing
    procedure, private :: assign_local_field
    procedure, private :: add_local_field
    procedure, private :: multiply
    generic :: operator(*) => multiply
    generic :: operator(+) => add_local_field
    generic :: assignment(=) => assign_local_field
    procedure :: output
  end type

  abstract interface
    pure function initial_condition(x) result(initial_values)
      import :: real64
      real(real64), intent(in) :: x
      real(real64) :: initial_values
    end function
  end interface

  real(real64), allocatable :: dx
  integer(int64), allocatable :: num_global_points,num_local_points
  integer(int64), parameter:: num_end_points=2_int64
  real(real64) :: boundary_vals(num_end_points)

contains

  function grid_spacing() result(delta_x)
    real(real64) :: delta_x
    call assert(allocated(dx),error_message("global_field%grid_spacing: dx not allocated"))
    delta_x = dx
  end function

  pure function state(this) result(local_values)
    class(global_field), intent(in) :: this
    real(real64), allocatable :: local_values(:)
    ! Requires
    if (this%user_defined()) then
      local_values = this%values
    end if
  end function

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

  subroutine set(this,initial_function,num_points)
    class(global_field), intent(inout) :: this
    integer, intent(in) :: num_points
    procedure(initial_condition), pointer :: initial_function
    integer(int64) :: num_intervals,i
    real(real64), parameter :: two_pi=2.*3.1415926535897932384626433832795028842_real64
    real(real64), allocatable :: local_grid(:)

    ! Requires
    call assert(mod(num_points,num_images())==0,error_message("global_field%set: num_points not evenly divisible by num_images()"))

    num_global_points=num_points
    num_local_points=num_points/num_images()
    num_intervals = num_global_points ! right-side boundary point is redundant and therefore not counted or stored
    dx=two_pi/real(num_intervals,real64)

    if (.not.allocated(this%values)) allocate(this%values(num_local_points)[*])
    local_grid = [((this_image()-1)*num_local_points+i-1,i=1,num_local_points)]*dx
    do concurrent(i=1:num_local_points)
      this%values(i) = initial_function(local_grid(i))
    end do
    call synchronize()

    ! Ensures
    call this%mark_as_defined
    call assert(allocated(dx),error_message("global_field%set: dx has not been allocated"))
    call assert(allocated(num_global_points),error_message("global_field%set: num_global_points has not been allocated"))
    call assert(allocated(num_local_points),error_message("global_field%set: num_local_points has not been allocated"))
  end subroutine

  subroutine assign_local_field(lhs,rhs)
    class(global_field), intent(inout) :: lhs
    class(local_field), intent(in) :: rhs
    real(real64), allocatable :: values(:)
    ! Requires
    if (.not.allocated(num_local_points)) error stop "global_field: no value established for memory allocation yet."
    if (.not.allocated(lhs%values)) allocate(lhs%values(num_local_points)[*])
    call assert(rhs%user_defined(),error_message("global_field%assign_local_field received uninitialized RHS."))
    lhs%values(:) = rhs%state()
    call synchronize()
    ! Ensures
    call lhs%mark_as_defined
  end subroutine

  pure function add_local_field(lhs,rhs) result(total)
    class(global_field), intent(in) :: lhs
    type(local_field), intent(in) :: rhs
    type(local_field) :: total
    ! Requires
    if (lhs%user_defined() .and. rhs%user_defined()) then
      total = lhs%values + rhs%state()
      call total%mark_as_defined
    end if
  end function

  pure function multiply(lhs,rhs) result(product_)
    class(global_field), intent(in) :: lhs,rhs
    type(local_field) :: product_
    ! Requires
    if (lhs%user_defined() .and. rhs%user_defined()) then
      product_= lhs%values * rhs%values
      call product_%mark_as_defined
    end if
  end function

 pure function x(this) result(this_x)
    class(global_field), intent(in) :: this
    type(local_field) :: this_x
    real(real64) :: local_this_x(num_local_points)
    integer(int64) :: i,left_neighbor,right_neighbor
    ! Requires
    if (this%user_defined() .and. allocated(dx) .and. allocated(num_local_points)) then
      associate(N=>num_local_points,me=>this_image())
        left_neighbor = merge(num_images(),me-1,me==1)
        local_this_x(1)=(this%values(2)-this%values(N)[left_neighbor])/(2._real64*dx)
        do concurrent(i=2:N-1)
          local_this_x(i)=(this%values(i+1)-this%values(i-1))/(2._real64*dx)
        end do
        right_neighbor = merge(1,me+1,me==num_images())
        local_this_x(N)=(this%values(1)[right_neighbor]-this%values(N-1))/(2._real64*dx)
      end associate
      this_x = local_this_x
      ! Ensures
      call this_x%mark_as_defined
    end if
  end function

 !pure function xx(this) result(this_xx)
  function xx(this) result(this_xx)
    class(global_field), intent(in) :: this
    type(local_field) :: this_xx
    real(real64) :: local_this_xx(num_local_points)
    integer(int64) :: i,left_neighbor,right_neighbor
    ! Requires
    if (this%user_defined() .and. allocated(dx) .and. allocated(num_local_points)) then
      associate(N=>num_local_points,me=>this_image())
        left_neighbor = merge(num_images(),me-1,me==1)
        local_this_xx(1)=(this%values(2)-2._real64*this%values(1)+this%values(N)[left_neighbor])/dx**2
        do concurrent(i=2:N-1)
          local_this_xx(i)=(this%values(i+1)-2._real64*this%values(i)+this%values(i-1))/dx**2
        end do
        right_neighbor = merge(1,me+1,me==num_images())
        local_this_xx(N)=(this%values(1)[right_neighbor]-2._real64*this%values(N)+this%values(N-1))/dx**2
      end associate
      this_xx = local_this_xx
      ! Ensures
      call this_xx%mark_as_defined
     !print *,"On image ",this_image(),", local_this_xx=",local_this_xx
     !stop
    end if
  end function

  subroutine output(this,unit,iotype,v_list,iostat,iomsg)
    class(global_field), intent(in) :: this
    integer, intent(in) :: unit ! Unit on which output happens (negative for internal file)
    character(*), intent(in) :: iotype ! Allowable values: ’LISTDIRECTED’,’NAMELIST’, or ’DT’
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    integer(int64) :: i
    ! Requires
    call assert(this%user_defined(),error_message("global_field%output received uninitialized object"))
    do i=1,size(this%values)
      write(unit,iostat=iostat) (this_image()-1)*size(this%values) + i, this%values(i)
    end do
  end subroutine

end module
