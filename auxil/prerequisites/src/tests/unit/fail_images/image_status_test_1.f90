! Check the status of all images. Error only, when one unexpectedly failed.

program test_image_status_1
  use iso_fortran_env , only : STAT_STOPPED_IMAGE
  implicit none
  integer :: i

  associate(np => num_images(), me => this_image())
    do i= 1, np
      if (image_status(i) /= 0) error stop "image_status(i) should not fail"
    end do

    sync all
    if (me == 1) print *,"Test passed."
  end associate

end program test_image_status_1

