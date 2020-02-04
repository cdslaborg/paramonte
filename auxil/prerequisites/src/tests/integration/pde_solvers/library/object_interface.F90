! MoFo library: object_interface
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

module object_interface
#include "compiler_capabilities.txt"
  implicit none
  private
  public :: object

  ! Define an abstract base class to ensure basic functionality expected to be provided by all concrete Morfeus classes.
  ! Each concrete class provides the functionality by extending this class and implementing its deferred binding(s).  This
  ! class resembles java's Object class in the sense that it is intended to be the ultimate ancester of every other class.
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
