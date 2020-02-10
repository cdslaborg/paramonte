! Unit test for register and allocated procedure.
!
! Test that matrix valued allocatable components in a derived typed coarray is
! registered correctly, delayed allocatable and deregisterable. The checks
! whether a component is allocated are done on this_image only.
!
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
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

program register_alloc_comp_3
  implicit none

  type dt
    integer, allocatable, dimension(:,:) :: m
  end type dt

  integer :: np = -2
  type(dt), allocatable :: obj[:]

  np = num_images()

  ! Allocate only the container. obj%i must not be allocated hereafter.
  allocate(obj[*])
  if (.not. allocated(obj)) error stop "Test failed. 'obj' not allocated."
  if (allocated(obj%m)) error stop "Test failed. 'obj%m' is allocated."
 
  ! Allocate the component.
  allocate(obj%m(3,6), source=this_image())

  ! Now both objects have to be allocated and obj%m(1:3,1:6) set to this_image()
  if (.not. allocated(obj)) error stop "Test failed. 'obj' not allocated."
  if (.not. allocated(obj%m)) error stop "Test failed. 'obj%m' not allocated."
  if (any( ubound(obj%m) /= [3, 6])) error stop "Test failed. ubound(obj%m) /= [3, 6]."
  if (any (obj%m(:,:) /= this_image())) error stop "Test failed. obj%m(:) /= this_image()."

  ! Deallocate the component.
  deallocate(obj%m)

  ! and test, that only the component is deallocated, but not the container.
  if (allocated(obj%m)) error stop "Test failed. 'obj%m' still allocated."
  if (.not. allocated(obj)) error stop "Test failed. 'obj' no longer allocated."

  ! Now deallocate the container, too.
  deallocate(obj)

  ! and check, that it worked.
  if (allocated(obj)) error stop "Test failed. 'obj' still allocated."

  ! Failing tests would make whole program error out.  Therefore it is save
  ! to print the pass message on image one, only.
  if (this_image() == 1) print *, "Test passed."
end program

