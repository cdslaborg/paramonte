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

! Unit tests for co_sum
program main
  use iso_fortran_env, only : error_unit
  use iso_c_binding, only : c_int,c_double
#ifdef USE_EXTENSIONS
  use opencoarrays
#endif
  implicit none
  logical :: co_sum_c_int_verified=.false.,co_sum_c_double_verified=.false.

#ifdef USE_EXTENSIONS
  if (this_image()==1) print *,"Using the extensions from the opencoarrays module."
#endif

  ! Verify collective sum of integer data by tallying image numbers
  c_int_co_sum: block
    integer(c_int) :: i,me
    me=this_image()
    sync all
    call co_sum(me)
    if (me==sum([(i,i=1,num_images())])) then
      co_sum_c_int_verified=.true.
    else
      write(error_unit,"(2(a,i2))") "co_broadcast with integer(c_int) argument fails with result (",me,") on image",this_image()
    end if
  end block c_int_co_sum

  ! Verify collective sum by calculuating pi
  c_double_co_sum: block
    real(c_double), parameter :: four=4._c_double,one=1._c_double,half=0.5_c_double
    real(c_double), save :: pi
    integer(c_int) :: i,points_per_image
    integer(c_int), parameter :: resolution=1024_c_int ! Number of points used in pi calculation
    integer(c_int) :: me
    me=this_image()
    ! Partition the calculation evenly across all images
    if (mod(resolution,num_images())/=0) then
      write(error_unit,"(a)") "number of images doesn't evenly divide into number of points"
      error stop
    end if
    points_per_image=resolution/num_images()
    associate(n=>resolution,my_first=>points_per_image*(me-1)+1,my_last=>points_per_image*me)
      pi = sum([ (four/(one+((i-half)/n)**2),i=my_first,my_last) ])/n
    end associate
    sync all
    ! Replace pi on each image with the sum of the pi contributions from all images
    call co_sum(pi)
    associate (pi_ref=>acos(-1._c_double),allowable_fractional_error=>0.000001_c_double)
      if (abs((pi-pi_ref)/pi_ref)<=allowable_fractional_error) then
        co_sum_c_double_verified=.true.
      else
        write(error_unit,*) "co_broadcast with real(c_double) argument fails with result (",pi,") result on image ",me
      end if
    end associate
  end block c_double_co_sum

  if (.not. all([co_sum_c_int_verified,co_sum_c_double_verified])) error stop
  ! Wait for every image to pass
  sync all
  if (this_image()==1) print *, "Test passed."
end program
