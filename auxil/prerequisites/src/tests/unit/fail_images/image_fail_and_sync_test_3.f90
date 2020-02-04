! Check that after an image has failed a sync images ends correctly.
! Image two is to fail, all others to continue.

program image_fail_and_sync_test_3
  use iso_fortran_env , only : STAT_FAILED_IMAGE
  implicit none
  integer :: i, syncAllStat

  associate(np => num_images(), me => this_image())
    if (np < 3) error stop "I need at least 3 images to function."
    do i= 1, np
      if (image_status(i) /= 0) error stop "image_status(i) should not fail"
    end do

    sync all
    if (me == 2) fail image
    sync all(STAT=syncAllStat)
    if (image_status(2) /= STAT_FAILED_IMAGE) error stop "Expected STAT_FAILED_IMAGE for image 2."

    if (me == 1) print *,"Test passed."
  end associate

end program image_fail_and_sync_test_3

