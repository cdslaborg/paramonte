! Copyright (c) 2012-2016, Sourcery, Inc.
! All rights reserved.
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

! Unit tests for co_broadcast and co_sum
program main
  use iso_fortran_env, only : error_unit
  use iso_c_binding, only : c_int,c_double,c_char
#ifdef USE_EXTENSIONS
  use opencoarrays
#endif
  implicit none
  integer(c_int) :: me
  ! Set test failure as the default result
  logical :: c_char_test_passes=.false.,c_int_test_passes=.false.,c_double_test_passes=.false.

  ! Store the executing image number
  me=this_image()

#ifdef USE_EXTENSIONS
  if (me==1) print *,"Using the extensions from the opencoarrays module."
#endif

  ! Verify broadcasting of character data from image 1
  c_char_co_broadcast: block
    character(kind=c_char,len=14), save :: string_received[*]
    character(kind=c_char,len=*), parameter :: string_sent=c_char_"Hello, world!"! Character test message
    if (me==1) string_received=string_sent
    sync all
    call co_broadcast(string_received,source_image=1)
    if (string_received/=string_sent) then
      write(error_unit,*) "Incorrect co_broadcast(",string_received,") on image",me
    else
      c_char_test_passes=.true.
    end if
  end block c_char_co_broadcast

  ! Verify broadcasting of integer data from image 1
  c_int_co_broadcast: block
    integer(c_int), save :: integer_received[*]
    integer(c_int), parameter :: integer_sent=12345_c_int ! Integer test message
    if (me==1) integer_received=integer_sent
    sync all
    call co_broadcast(integer_received,source_image=1)
    if (integer_received/=integer_sent) then
      write(error_unit,*) "Incorrect co_broadcast(",integer_received,") on image",me
    else
      c_int_test_passes=.true.
    end if
  end block c_int_co_broadcast

  ! Verify broadcasting of real data from image 1
  c_double_co_broadcast: block
    real(c_double), save :: real_received[*]
    real(c_double), parameter :: real_sent=2.7182818459045_c_double ! Real test message
    if (me==1) real_received=real_sent
    sync all
    call co_broadcast(real_received,source_image=1)
    if (real_received/=real_sent) then
      write(error_unit,*) "Incorrect co_broadcast(",real_received,") on image",me
    else
      c_double_test_passes=.true.
    end if
  end block c_double_co_broadcast


  if (.not.all([c_char_test_passes,c_int_test_passes,c_double_test_passes])) error stop
  ! Wait for everyone to pass the tests
  sync all
  if (me==1) print *, "Test passed."
end program
