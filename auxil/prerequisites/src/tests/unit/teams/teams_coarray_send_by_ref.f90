program teams_coarray_get_by_ref
  use, intrinsic :: iso_fortran_env, only: team_type
  implicit none
  type(team_type) :: team
  type :: allocatable_array_t
    integer, allocatable :: A(:)
  end type
  type(allocatable_array_t) :: R[*]
  integer :: i, my_team, initial_team_this_image

  ! handle odd or even number of images
  allocate(R%A(num_images()/2+mod(num_images(),2)*mod(this_image(),2)), source=0)

  initial_team_this_image = this_image()
  my_team = mod(this_image()-1,2)+1

  form team (my_team, team)

  change team (team)
    do i = 1, num_images()
       R[i]%A(this_image()) = initial_team_this_image
    end do
  end team

  if (any(R%A /= [(i, i=my_team, num_images(), 2)])) error stop 'Test failed.'

  sync all

  if (this_image() == 1) write(*,*) 'Test passed.'
end program teams_coarray_get_by_ref
