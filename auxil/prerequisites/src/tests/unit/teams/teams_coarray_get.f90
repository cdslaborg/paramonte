program teams_coarray_get
  use, intrinsic :: iso_fortran_env, only: team_type
  implicit none
  type(team_type) :: team
  integer, allocatable :: L(:)
  integer :: i, my_team, R[*]

  ! handle odd or even number of images
  allocate(L(num_images()/2+mod(num_images(),2)*mod(this_image(),2)))

  R = this_image()
  my_team = mod(this_image()-1,2)+1

  form team (my_team, team)

  change team (team)
    do i = 1, num_images()
      L(i) = R[i]
    end do
  end team

  if (any(L /= [(i, i=my_team, num_images(), 2)])) error stop 'Test failed.'

  sync all

  if (this_image() == 1) write(*,*) 'Test passed.'
end program teams_coarray_get
