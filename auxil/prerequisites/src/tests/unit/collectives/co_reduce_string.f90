! Implement a usecase for co_reduce on char arrays.
!
! Copyright (c) 2012-2017, Sourcery, Inc.
! All rights reserved.
!
! Unit tests for co_min: verify parallel, collective minimum
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!     * Neither the name of Sourcery, Inc., nor the
!       names of any other contributors may be used to endorse or promote products
!       derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL SOURCERY, INC., BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

program co_reduce_strings
 
  implicit none

  integer, parameter :: numstrings = 10, strlen = 6
  character(len=strlen), dimension(:), allocatable :: strarr[:]
  character(len=strlen) :: expect
  integer :: i

  ! Construct the strings by postfixing foo by a number.
  associate (me => this_image())
    allocate(strarr(numstrings)[*])
    do i = 1, numstrings
      write(strarr(i), "('foo'I02)") i * me
    end do
    ! Collectively reduce the maximum string.
    call co_reduce(strarr, strmax)
  end associate

  ! No sync should be needed here, because the collective (reduce_all)
  ! implicitly synchronizes.
  associate (np => num_images())
    do i = 1, np
      write (expect, "('foo'I02)") i * np
      if (strarr(i) /= expect) then
        ! On errror print what we got and what we expected.
        print *, "Got: ", strarr(i), ", expected: ", expect
        error stop "Didn't get expected string."
      end if
    end do
  end associate
  sync all
  if (this_image() == 1) print *, "Test passed."
contains

  !! Compare two strings and return the maximum one. In a co_reduce no deferred-
  !! length strings are allowed, therefore fixed length had to be used.
  !! For identical strings the LHS is returned.
  pure function strmax(lhs, rhs) result(maxstr) bind(C,name="strmax")
    character(len=strlen), intent(in) :: lhs,rhs
    character(len=strlen) :: maxstr

    if (lhs > rhs) then
      maxstr = lhs 
    else 
      maxstr = rhs
    end if
  end function

end program co_reduce_strings

