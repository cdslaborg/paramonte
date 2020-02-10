! Check that after an image has failed a get to other images is still possible.
! Image two is to fail, all others to continue.

program image_fail_and_get_test_1
  use iso_fortran_env , only : STAT_FAILED_IMAGE
  implicit none
  integer :: i, stat
  integer, save :: share[*]

  associate(np => num_images(), me => this_image())
    if (np < 3) error stop "I need at least 3 images to function."
    
    share = 37

    sync all
    if (me == 2) fail image
    sync all (STAT=stat)

    print *, "Checking shared value."
    do i= 1, np
      if (i /= 2 .AND. i /= me) then
        if (share[i, STAT=stat] /= 37) error stop "Expected to get value from images alive."
        print *, me, "Stat of #", i, " is:", stat
      end if
    end do

    sync all(STAT=stat)
    if (me == 1) print *,"Test passed."
  end associate

end program image_fail_and_get_test_1

