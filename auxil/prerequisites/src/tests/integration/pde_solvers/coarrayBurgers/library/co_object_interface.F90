module co_object_interface
  implicit none
  private
  public :: co_object

  ! Define an abstract base class to ensure basic functionality expected to be provided by all concrete Morfeus classes.
  ! Each concrete class provides the functionality by extending this class and implementing its deferred binding(s).  This
  ! class resembles java's Object class in the sense that it is intended to be the ultimate ancester of every other class.
  type, abstract :: co_object
    private
    logical :: defined=.false. ! Mark all co_objects as not-yet user-defined by default
    real, allocatable :: dummy_to_facilitate_extension[:]
  contains
    procedure :: mark_as_defined
    procedure :: user_defined
    procedure(formatted_output_interface), deferred :: output
   !generic :: write(unformatted) => output  ! Derived-type I/O not yet supported by most compilers
  end type

  ! Require child classes to write an "output" procedure that prints to the passed file unit
  abstract interface
    subroutine formatted_output_interface(this,unit,iotype,v_list,iostat,iomsg)
      import co_object
      class(co_object), intent(in) :: this
      integer, intent(in) :: unit ! Unit on which output happens (negative for internal file)
      character(*), intent(in) :: iotype ! Allowable values: ’LISTDIRECTED’,’NAMELIST’, or ’DT’
      integer, intent(in) :: v_list(:)
      integer, intent(out) :: iostat
      character(*), intent(inout) :: iomsg
    end subroutine
  end interface

contains

  ! Mark the co_object as user-defined
  pure subroutine mark_as_defined(this)
    class(co_object), intent(inout) :: this
    this%defined=.true.
  end subroutine

  ! Return a boolean result indicating whether this co_object has been initialized since its declaration
  logical pure function user_defined(this)
    class(co_object), intent(in) :: this
    user_defined = this%defined
  end function

end module
