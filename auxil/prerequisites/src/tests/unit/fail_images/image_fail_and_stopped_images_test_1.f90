! Check that letting images exit the stopped_images function returns the
! correct indices.
! Image two is to stop, all others to continue.

program image_fail_and_stopped_images_test_1
  implicit none
  integer :: i, stat
  integer, allocatable :: simages(:)

  associate(np => num_images(), me => this_image())
    if (np < 3) error stop "I need at least 3 images to function."

    ! Need a sync here to make sure all images are started and to prevent image fail
    ! is not already detected in above image_status(i).
    sync all
    if (me == 2) stop 0
    sync all (STAT=stat)

    simages = stopped_images()
    if (size(simages) /= 1) error stop "stopped_images()'s size should be one."
    if (simages(1) /= 2) then
            print *, me, "stopped image: ", simages(1)
            error stop "The second image should have stopped."
    end if

    if (me == 1) print *,"Test passed."
  end associate

end program image_fail_and_stopped_images_test_1

