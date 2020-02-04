module local_field_module
   use iso_fortran_env, only : real64,int64
   use ForTrilinos_assertion_utility, only : assert,error_message
   use object_interface, only : object
   implicit none
   private
   public :: local_field

   type, extends(object) :: local_field
     private
     real(real64), allocatable :: values(:)
   contains
     procedure :: state
     procedure, private, pass(rhs) :: multiply
     procedure, private :: subtract
     procedure, private :: assign_array
     generic :: operator(-)=>subtract
     generic :: operator(*)=>multiply
     generic :: assignment(=)=>assign_array
     procedure :: output
   end type

contains

  pure subroutine assign_array(lhs,rhs)
    class(local_field), intent(inout) :: lhs
    real(real64), intent(in) :: rhs(:)
    lhs%values = rhs
    ! Ensures
    call lhs%mark_as_defined
  end subroutine

  pure function subtract(lhs,rhs) result(difference)
    class(local_field), intent(in) :: lhs,rhs
    type(local_field) :: difference
    !Requires
    if (lhs%user_defined() .and. rhs%user_defined()) then
      difference%values = lhs%values - rhs%values
      ! Ensures
      call difference%mark_as_defined
    end if
  end function

  pure function multiply(lhs,rhs) result(product_)
    class(local_field), intent(in) :: rhs
    type(local_field) :: product_
    real(real64), intent(in) :: lhs
    if (rhs%user_defined()) then
      product_%values = lhs*rhs%values
      ! Ensures
      call product_%mark_as_defined
    end if
  end function

  pure function state(this) result(this_values)
    class(local_field), intent(in) :: this
    real(real64), allocatable :: this_values(:)
    this_values = this%values
  end function

  subroutine output(this,unit,iotype,v_list,iostat,iomsg)
    class(local_field), intent(in) :: this
    integer, intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg
    integer(int64) :: i
    ! Requires
    call assert(this%user_defined(),error_message("local_field%output received uninitialized object"))
    do i=1,size(this%values)
      write(unit,iostat=iostat) (this_image()-1)*size(this%values) + i, this%values(i)
    end do
  end subroutine

end module
