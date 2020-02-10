! MoFo library: parse_command_line
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

module parse_command_line
  ! Utility for returning key-value pairs passed at the command line in the format
  ! command key1=value1  key2=value2 ...
  implicit none
  private
  public :: get_keyword_values

contains

  subroutine assert(test_passes)
    logical, intent(in) :: test_passes
    if (.not. test_passes) stop "assertion failed"
  end subroutine

  ! Return the command-line values associated with the passed keys.
  ! Arguments and result variables:
  ! keys = array of keywords
  ! default_values = values returned for like-positioned keyword if no command-line value specified
  ! actual_values = values returned for like-positioned keys
  ! Each elment in the keys and default_values arrays must be padded with trailing blanks if necessary to give the elements
  ! a uniform length.  The actual_values elements are similarly padded if necessary to give them a uniform length.
  function get_keyword_values(keys,default_values) result(actual_values)
    character(len=*), dimension(:), intent(in) :: keys
    character(len=*), dimension(:), intent(in) :: default_values
    character(len=:), dimension(:), allocatable :: actual_values
    character(len=:), allocatable :: key_value_pair,trimmed_argument,trimmed_value
    character(len=1), parameter :: divider="="
    integer divider_position,error_flag,i,j

    ! Requires
    call assert(size(keys)==size(default_values))

    actual_values=default_values
    allocate(key_value_pair,source=repeat(" ",ncopies=max_argument_length()) )
    ! Read the text of the arguments passed on the command line
    do i=1,command_argument_count()
      call get_command_argument(i,key_value_pair,status=error_flag)
      call check(error_flag)
      divider_position = scan(key_value_pair,divider)
      if (divider_position==0) stop "Invalid argument format (expected: 'argument=value')."
      trimmed_argument = trim(key_value_pair(1:divider_position-1))
      trimmed_value = trim(key_value_pair(divider_position+1:))
      if (len(trimmed_value)==0) stop "Invalid value format (expected: 'argument=value')."
      do j=1,size(keys)
        if (trim(keys(j))==trimmed_argument) actual_values(j)=trimmed_value
      end do
    end do

    ! Ensures
    call assert(size(actual_values)==size(keys))
  contains
    function max_argument_length()
      integer max_argument_length,n,length_of_argument_n
      max_argument_length=0
      do n=1,command_argument_count()
        call get_command_argument(n,key_value_pair,status=error_flag,length=length_of_argument_n)
        call check(error_flag)
        if (length_of_argument_n>max_argument_length) max_argument_length = length_of_argument_n
      end do
    end function
    subroutine check(flag)
      integer, intent(in) :: flag
      select case(flag)
        case(-1)
          ! this should never occur because key_value_pair is dynamically sized to match the length of the longest argument
          print *,"main: argument ",i,"exceeds maximum length of ",max_argument_length()
          stop
        case(1:)
          print *,"main: error in reading the argument name (status=",flag,")"
          stop
        case(0)
          ! argument_text read successfully
        case default
          stop "main: invalid status (compiler error)"
      end select
    end subroutine
  end function
end module
