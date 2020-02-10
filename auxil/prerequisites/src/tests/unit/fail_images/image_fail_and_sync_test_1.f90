! Check that after an image has failed a sync all ends correctly.
! Image two is to fail, all others to continue.

program image_fail_and_sync_test_1
  use iso_fortran_env , only : STAT_FAILED_IMAGE
  implicit none
  integer :: i, syncAllStat

  associate(np => num_images(), me => this_image())
    if (np < 3) error stop "I need at least 3 images to function."
    do i= 1, np
      if (image_status(i) /= 0) error stop "image_status(i) should not fail"
    end do

    ! Need a sync here to make sure all images are started and to prevent image fail
    ! is not already detected in above image_status(i).
    sync all
    if (me == 2) fail image
    sync all (STAT=syncAllStat)
    ! Check that all images returning from the sync report the failure of an image
    if (syncAllStat /= STAT_FAILED_IMAGE) error stop "Expected sync all (STAT == STAT_FAILED_IMAGE)."

    if (me == 1) print *,"Test passed."
  end associate

end program image_fail_and_sync_test_1

