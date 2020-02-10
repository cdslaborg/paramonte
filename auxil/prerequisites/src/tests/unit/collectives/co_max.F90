! Copyright (c) 2012-2016, Sourcery, Inc.
! All rights reserved.
!
! Unit tests for co_max: verify parallel, collective maximum
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

program main
  use iso_fortran_env, only : error_unit
  use iso_c_binding, only : c_int,c_double
#ifdef USE_EXTENSIONS
  use opencoarrays
#endif
  implicit none

#ifdef USE_EXTENSIONS
  if (this_image()==1) print *,"Using the extensions from the opencoarrays module."
#endif

  ! Verify that 1 is the lowest image number
  c_int_co_min: block
    integer(c_int) :: me
    me=this_image()
    sync all
    call co_max(me)
    if (me/=num_images()) then
      write(error_unit,"(2(a,i2))") "Wrong result (",me,") on image",this_image()
      error stop
    end if
    ! Wait for all images to pass the test
    sync all
    if (me==1) print *,"Correct integer(c_int) co_min"
  end block c_int_co_min

  ! Verify that the maximum real conversion of an image number is real(num_images())
  c_double_co_max: block
    real(c_double) :: me
    me=real(this_image(),c_double)
    sync all
    call co_max(me)
    if (me/=real(num_images(),c_double)) then
      write(error_unit,"(2(a,i2))") "Wrong result (",me,") on image",this_image()
      error stop
    end if
    ! Wait for all images to pass the test
    sync all
    if (this_image()==1) print *,"Correct real(c_double) co_max"
  end block c_double_co_max

  if (this_image()==1) print *, "Test passed."
end program
