program teams_coarray_send
  use, intrinsic :: iso_fortran_env, only: team_type
  implicit none
  type(team_type) :: team
  integer, allocatable :: R(:)[:]
  integer :: extent, i, my_team, initial_team_this_image, odd

  ! if odd number of images, even images have R(extent) == 0
  extent = num_images()/2+mod(num_images(),2)
  allocate(R(extent)[*], source=0) 

  initial_team_this_image = this_image()
  my_team = mod(this_image()-1,2)+1

  form team (my_team, team)

  change team (team)
    do i = 1, num_images()
      R(this_image())[i] = initial_team_this_image
    end do
  end team

  if (any(R /= [(mod(i, num_images()+1), i=my_team, 2*extent, 2)])) error stop 'Test failed.'

  sync all

  if (this_image() == 1) write(*,*) 'Test passed.'
end program teams_coarray_send
