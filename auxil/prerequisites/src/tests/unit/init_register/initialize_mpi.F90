! Unit test for initializion of MPI by LIBCAF_MPI.
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
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
program initialize_mpi
#ifdef MPI_WORKING_MODULE
  use mpi, only : MPI_COMM_SIZE,MPI_COMM_WORLD
  implicit none
#else
  implicit none
  include 'mpif.h'
  interface
     subroutine MPI_COMM_SIZE(mpi_comm,nranks,ierr)
       integer, intent(in)  :: mpi_comm
       integer, intent(out) :: nranks, ierr
     end subroutine
  end interface
#endif

  ! Set invalid default image number and number of ranks
  integer :: me=-1,np=-1,ierr

  ! Get image number
  me = this_image()

  ! Get number of ranks (np)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

  ! Everybody verifies that they have a valid image number and rank
  if(me < 1 .or. np < 1 .or. me > np) error stop "Test failed."

  sync all

  ! Image 1 reports test success
  if(me==1) print *,"Test passed."
end program
