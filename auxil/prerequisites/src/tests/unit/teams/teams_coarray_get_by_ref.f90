program teams_coarray_get_by_ref
  use, intrinsic :: iso_fortran_env, only: team_type
  implicit none
  type(team_type) :: team
  type :: allocatable_array_t
    integer, allocatable :: A(:)
  end type
  type(allocatable_array_t) :: R[*]
  integer, allocatable :: L(:)
  integer :: i, my_team

  ! handle odd or even number of images
  allocate(L(num_images()/2+mod(num_images(),2)*mod(this_image(),2)))

  my_team = mod(this_image()-1,2)+1

  form team (my_team, team)

  ! size(R%A) == this_image(team)
  allocate(R%A((this_image()+1)/2), source=0)

  R%A(ubound(R%A,1)) = this_image()

  change team (team)
    do i = 1, num_images()
      L(i) = R[i]%A(i)
    end do
  end team

  if (any(L /= [(i, i=my_team, num_images(), 2)])) error stop 'Test failed.'

  sync all

  if (this_image() == 1) write(*,*) 'Test passed.'
end program teams_coarray_get_by_ref
