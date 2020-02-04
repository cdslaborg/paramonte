module object_interface
#include "compiler_capabilities.txt"
  implicit none
  private
  public :: object

  ! Define an abstract parent type to ensure basic functionality expected to be provided by all non-abstract types.
  ! Each non-abstract type provides the functionality by extending this type and implementing its deferred binding(s).  This
  ! type resembles java's Object class in the sense that it is intended to be the ultimate ancester of every other type.
  type, abstract :: object
    private
    logical :: defined=.false. ! Mark all objects as not-yet user-defined by default
  contains
    procedure :: mark_as_defined
    procedure :: user_defined
    procedure(output_interface), deferred :: output
#ifndef COMPILER_LACKS_DERIVED_TYPE_IO
    generic :: write(formatted) => output  ! Derived-type I/O
#endif /* COMPILER_LACKS_DERIVED_TYPE_IO */
  end type

  ! Require child classes to write an "output" procedure that prints to the passed file unit
  abstract interface

    subroutine output_interface(this,unit,iotype,v_list,iostat,iomsg)
      import object
      class(object), intent(in) :: this
      integer, intent(in) :: unit
      character(len=*), intent(in) :: iotype
      integer, intent(in) :: v_list(:)
      integer, intent(out) :: iostat
      character(len=*), intent(inout) :: iomsg
    end subroutine

  end interface

contains

  ! Mark the object as user-defined
  pure subroutine mark_as_defined(this)
    class(object), intent(inout) :: this
    this%defined=.true.
  end subroutine

  ! Return a boolean result indicating whether this object has been initialized since its declaration
  logical pure function user_defined(this)
    class(object), intent(in) :: this
    user_defined = this%defined
  end function

end module
