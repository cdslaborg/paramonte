! BSD 3-Clause License
!
! Copyright (c) 2018, Sourcery Institute
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
  !! summary: Test get_commiunicator function, an OpenCoarrays-specific language extension
  use opencoarrays, only : get_communicator
  use oc_assertions_interface, only : assert

  implicit none

#ifndef MPI_WORKING_MODULE
  include 'mpif.h'
#endif

  call mpi_matches_caf(get_communicator())
    !! verify # ranks = # images and image number = rank + 1

  block
    use iso_fortran_env, only : team_type
    use opencoarrays, only : get_communicator, team_number !! TODO: remove team_number once gfortran supports it

    type(team_type) :: league
    integer, parameter :: num_teams=2
      !! number of child teams to form from the parent initial team

    associate(initial_image=>this_image(), initial_num_images=>num_images(), chosen_team=>destination_team(this_image(),num_teams))

      form team(chosen_team,league)
        !! map images to num_teams teams

      change team(league)
        !! join my destination team

        call mpi_matches_caf(get_communicator())
          !! verify new # ranks = new # images and new image number = new rank + 1

        associate(my_team=>team_number())
          call assert(my_team==chosen_team,"assigned team matches chosen team")
          associate(new_num_images=>initial_num_images/num_teams+merge(1,0,my_team<=mod(initial_num_images,num_teams)))
           call assert(num_images()==new_num_images,"block distribution of images")
          end associate
        end associate

      end team

      call assert( initial_image==this_image(),"correctly remapped to original image number")
      call assert( initial_num_images==num_images(),"correctly remapped to original number of images")

    end associate

  end block

  sync all
  if (this_image()==1) print *,"Test passed."

contains

   pure function destination_team(image,numTeams) result(team)
     integer, intent(in) ::image, numTeams
     integer ::team
     team = mod(image+1,numTeams)+1
   end function

subroutine mpi_matches_caf(comm)
#ifdef MPI_WORKING_MODULE
    use mpi
#endif
    use iso_c_binding, only : c_int
    integer(c_int), intent(in) :: comm
      !! MPI communicator
    integer(c_int) :: isize,ierror,irank

    call MPI_COMM_SIZE(comm, isize, ierror)
    call assert( ierror==0 , "successful call MPI_COMM_SIZE" )
    call assert( isize==num_images(), "num MPI ranks = num CAF images " )

    call MPI_COMM_RANK(comm, irank, ierror)
    call assert( ierror==0 , "successful call MPI_COMM_RANK" )
    call assert( irank==this_image()-1 , "correct rank/image-number correspondence" )

  end subroutine

end program
