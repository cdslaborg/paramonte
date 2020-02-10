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

program main
  !! author: Damian Rouson
  !! date: 2017-08-01
  !! category: regression
  !! Test whether assigning an allocatable coarray component to an allocatable
  !! coarray component works.
  !! OpenCoarrays issue #422
  use iso_fortran_env, only : error_unit
  implicit none

  type foo
    !! A derived type is required to demonstrate issue #422
    integer, allocatable :: bar(:)[:]
  end type
  type(foo) :: foobar

  enum, bind(C)
    enumerator :: recipient=1, provider
    !! provider=2
  end enum

  integer, parameter :: message(*)=[3,4], required_images=2
    !! Data for provider image to put on recipient image

  associate( me=>this_image(), N_images=>num_images() )

    verify_num_images: if (N_images<required_images) then
      write(error_unit,*)  "issue-422-send.f90 requires at least ",required_images," images"
      error stop
    end if verify_num_images

    allocate(foobar%bar(size(message))[*],source=me)
      !! Assign each image's identifier to each element of foobar%bar
#ifndef GCC_GE_7
    sync all
      !! Fix to #243 has not been backported to GCC < 7. Implicit sync
      !! happens after allocation but before sourced assignment
#endif

    get_put_component: if (me==provider) then
      foobar%bar(:)[recipient] =  message
        !! Get bar from provider image and put bar on recipient image
      sync images(recipient)
        !! Signal recipient that get and put have completed
    end if get_put_component

    wait_and_verify: if (me==recipient) then
      sync images(provider)
        !! Wait for signal from provider image

      verify_result: if (any(foobar%bar/=message)) then
        write(error_unit,*) "Recipient image ",recipient," received ",foobar%bar," but expected ",provider
        error stop
      end if verify_result

      print *,"Test passed."
        !! Report success

    end if wait_and_verify

  end associate

end program
