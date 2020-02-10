! async_comp_alloc.f90 

! Unit test for register procedure and remote allocated test:
! Test that scalar allocatable component in a derived typed coarray is
! registered correctly, its allocation is deferred until the program
! allocates it specifically, and it is deregisterable. The component 
! allocation checks are done on remote images.
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

program async_comp_alloc
  implicit none

  type dt
    integer, allocatable :: i
  end type dt
  type(dt), allocatable :: obj[:]
    !! Container object for testing component allocation status

  logical, parameter :: verbose=.true. 
    !! Toggle verbose output for debugging purposes

  integer :: allocation_status
  integer, parameter :: success=0
    !! Successful allocation status value


  associate(me=>this_image(),np=>num_images())
    !! Make me & np locally immutable by associating them with function results

    call assert( np > 1, "num_images()>1")
      !! Ensure at least two images to simplify the test 

    allocate(obj[*],stat=allocation_status)
      !! Allocate only the container. obj%i must not be allocated hereafter.

    call assert(allocation_status == success, "allocated(obj)")
    call assert(.not. allocated(obj%i)      , ".not. allocated(obj%i)")

    block 
      integer :: allocating_image, test_image
      character(len=20) :: image_number

      loop_over_all_image_numbers: do allocating_image = 1, np
        !! Check that all allocations have been performed up to allocating_image and
        !! that no other allocations have been performed
  
        if (verbose) print *, me, "/", np, ": allocating_image=", allocating_image
  
        if (allocating_image /= me) then
          sync all
           !! Order allocations sequentially by image number: arrival of image 'allocating_image' 
           !! at the sync all in the alternate branch launches the following test loop.
          test_allocation_status: do test_image = 1, np
            if (verbose) print *, me, "/", np, ": Checking", test_image, " for allocation status."
            write(image_number,*) test_image
            if (test_image <= allocating_image) then
              call assert(allocated(obj[test_image]%i), " allocated(obj["//image_number//"]%i)" )
                !! Enforce that image numbers up to mine have allocated their components
            else
              call assert(.not. allocated(obj[test_image]%i), ".not. allocated(obj["//image_number//"]%i) on image ")
                !! Enforce that image numbers higher than 'image' should not have allocated their components yet
            end if
          end do test_allocation_status 
        else
          if (verbose) print *, me, "/", np, ": allocating..."
          allocate(obj%i, source = me)
            !! TODO: should we also ensure that implicit (re)allocation upon assignment works here?
          if (verbose) print *, me, "/", np, ": allocated"
          ! Enforce that object has been allocated and obj%i is this_image():
          call assert( allocated(obj)  , "allocated(obj)")
          call assert( allocated(obj%i), "allocated(obj%i)")
          call assert( obj%i == me     , "obj%i == this_image()")
          sync all
        end if
        sync all
      end do loop_over_all_image_numbers
    end block
    if (me == 1) print *, "Test passed."
  end associate
contains
  subroutine assert(assertion,description,stat)
    logical, intent(in) :: assertion 
    character(len=*), intent(in) :: description
    integer, intent(out), optional:: stat
    integer, parameter :: failure=1
    if (assertion) then
       if (present(stat)) stat=success
    else
      if (present(stat)) then
        stat=failure
      else
        error stop "Assertion "// description //" failed."
      end if
    end if
  end subroutine
end program

