! test1caf test
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
!

program test1caf
  implicit none
  integer, parameter :: num_local_elems=3,a_initial=1,b_initial=2
  integer :: a(num_local_elems)[*]=a_initial,b(num_local_elems)[*]=b_initial
  integer :: i,me,np,left,right

  me = this_image()
  np = num_images()

  left  = merge(np,me-1,me==1)
  right = merge(1,me+1,me==np)

  if (mod(me,2).eq.0) then
     a(:)[right] = a(:)[right]+me
  else
     b(:)[left] = b(:)[left]+me
  end if

  if(me==1) then
     write(*,*) me, a, b
  else
     sync images(me-1)
     write(*,*) me, a, b
  end if

  if(me < np) sync images(me+1)

  sync all

  if (mod(me,2).eq.0) then
    if ( any(a(:)[right]/=a_initial+me)) error stop "Test failed."
  else
    if ( any(b(:)[left]/=b_initial+me)) error stop "Test failed."
  end if

  sync all

  if (me==1) print *,"Test passed."

end program test1caf
