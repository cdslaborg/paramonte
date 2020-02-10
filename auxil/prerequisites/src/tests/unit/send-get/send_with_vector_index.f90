program send_with_vector_index

  implicit none

  integer, parameter :: indexes(4) = [ 2, 4, 1, 3]
  integer, dimension(4) :: dst[*]

  associate (me => this_image(), np => num_images())
    if (np < 2) error stop "Need at least two images."

    dst = -1
    sync all

    if (me == 2) dst(indexes)[1] = [ 99, 999, 9999, 99999]

    sync all

    if (me == 1) then
      print *, "me=", me, ":", dst
      if (any(dst /= [9999, 99, 99999, 999])) error stop "Test failed."
    else if (me == 2) then
      dst = [ 2, 33, 222, 333]
      dst(indexes)[2] = dst
      print *, "me=", me, ":", dst
      if (any(dst /= [222, 2, 333, 33])) error stop "Test failed."
    end if

    sync all
    if (me == 1) print *, 'Test passed.'
  end associate
end program send_with_vector_index

! vim:ts=2:sts=2:sw=2:
