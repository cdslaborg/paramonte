! BSD 3-Clause License
!
! Copyright (c) 2016, Sourcery Institute
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! * Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!
! * Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
!
! * Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! Comments preceded by "!!" are formatted for the FORD docoumentation generator
program allocatable_p2p_event_post
  !! author: Andre Vehreschild
  !! date: 2017-08-05
  !! category: unit-test
  !! Basic events test testing receipt of event post from one image to another
  use iso_fortran_env, only: event_type
  implicit none

  type(event_type), allocatable :: snd_copied(:)[:]

  if (num_images() < 4) error stop "num_images() >= 4 required for even_post_1 test"
  associate(me => this_image(), np => num_images())
    allocate(snd_copied(np)[*])
    if (me == 2) print*,'I am  image 2, I am posting to 4'
    if (me == 2) event post(snd_copied(2)[4])
    if (me == 2) print*,' I am image 2, I have posted to 4'
    if (me == 4) then
      event wait(snd_copied(2))
      sync all ! sync not required, but *may* expose cleanup issues/segfaults etc.
      print *, 'Test passed.'
    end if
    if (me /= 4) then
      sync all
      print *, 'I am', me, 'and image 4 told me it received the event'
    end if
    if (allocated(snd_copied)) deallocate(snd_copied)
  end associate
end program
! vim:ts=2:sts=2:sw=2:
