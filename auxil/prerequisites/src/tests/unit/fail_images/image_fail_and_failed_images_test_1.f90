! Check that after an image has failed the failed_images function returns the
! correct indices.
! Image two is to fail, all others to continue.

program image_fail_and_failed_images_test_1
  use iso_fortran_env , only : STAT_FAILED_IMAGE
  implicit none
  integer :: i, stat
  integer, allocatable :: fimages(:)

  associate(np => num_images(), me => this_image())
    if (np < 3) error stop "I need at least 3 images to function."
    do i= 1, np
      if (image_status(i) /= 0) error stop "image_status(i) should not fail"
    end do

    ! Need a sync here to make sure all images are started and to prevent image fail
    ! is not already detected in above image_status(i).
    sync all
    if (me == 2) fail image
    sync all (STAT=stat)

    fimages = failed_images()
    if (size(fimages) /= 1) error stop "failed_images()'s size should be one."
    if (fimages(1) /= 2) error stop "The second image should have failed."

    if (me == 1) print *,"Test passed."
  end associate

end program image_fail_and_failed_images_test_1

