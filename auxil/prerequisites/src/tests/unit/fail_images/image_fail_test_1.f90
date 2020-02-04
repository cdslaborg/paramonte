! Check that failing an image works.
! Image two is to fail, all others to continue.

program image_fail_test_1
  use iso_fortran_env , only : STAT_FAILED_IMAGE
  implicit none
  integer :: i, syncAllStat

  associate(np => num_images(), me => this_image())
    if (np < 3) error stop "I need at least 3 images to function."

    if (me == 2) fail image

    if (me == 2) print *, "Test failed."

    sync all (STAT=syncAllStat)
    if (me == 1) print *, "Test passed."
  end associate

end program image_fail_test_1
