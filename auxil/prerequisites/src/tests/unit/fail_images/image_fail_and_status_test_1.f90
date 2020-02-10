! Check that the status of a failed image is retrieved correctly.
! Image two is to fail, all others to continue.

program image_fail_and_status_test_1
  use iso_fortran_env , only : STAT_FAILED_IMAGE, STAT_STOPPED_IMAGE
  implicit none
  integer :: i, stat

  associate(np => num_images(), me => this_image())
    if (np < 3) error stop "I need at least 3 images to function."
    do i= 1, np
      if (image_status(i) /= 0) error stop "image_status(i) should not fail"
    end do

    ! Need to sync here or above image_status might catch a fail already.
    sync all
    if (me == 2) fail image
    sync all (STAT=stat)
    ! Check that all images returning from the sync report the failure of an image
    print *,"sync all (STAT=", stat, ")"
    if (stat /= STAT_FAILED_IMAGE) error stop "Expected sync all (STAT == STAT_FAILED_IMAGE)."

    do i= 1, np
      stat = image_status(i)
      if (i /= 2 .AND. stat /= 0 .AND. stat /= STAT_STOPPED_IMAGE) error stop "image_status(i) should not fail"
      if (i == 2 .AND. stat /= STAT_FAILED_IMAGE) error stop "image_status(2) should report fail"
    end do

    if (me == 1) print *,"Test passed."
  end associate

end program image_fail_and_status_test_1

