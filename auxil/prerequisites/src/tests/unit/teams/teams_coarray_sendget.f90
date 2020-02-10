program teams_coarray_sendget
  use, intrinsic :: iso_fortran_env, only: team_type
  implicit none
  type(team_type) :: team
  integer, allocatable :: R_send(:,:)[:]
  integer :: extent, i, j, my_team, team_num_images, R_get[*]

  ! if there are an odd number of images, then even images have R(:,extent) == 0
  extent = num_images()/2+mod(num_images(),2)
  allocate(R_send(extent, extent)[*], source=0)

  my_team = mod(this_image()-1,2)+1

  form team (my_team, team)

  R_get = this_image()

  change team (team)
    team_num_images = num_images()
    do concurrent (i = 1:num_images(), j = 1:num_images())
      R_send(this_image(),j)[i] = R_get[j]
    end do
  end team

  if (any(R_send /= reshape([((merge(i,0,i<=num_images() .and. j <= team_num_images), &
                              i=my_team,2*extent,2),j=1,extent)], &
                            shape=[extent,extent], order=[2,1]))) error stop 'Test failed.'

  sync all

  if (this_image() == 1) write(*,*) 'Test passed.'
end program teams_coarray_sendget
