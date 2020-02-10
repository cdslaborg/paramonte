program sync_image_ring_abort_on_stopped_image
  !! `SYNC IMAGES([this_image - 1, this_image + 1])` with
  !! `STAT=STAT_STOPPED_IMAGE` specifier on a periodic ring.  The test
  !! checks that syncing in a ring with a stopped image still
  !! terminates all images. All images other than image 1 participate
  !! in the `sync images()` call

  use, intrinsic:: iso_fortran_env
  implicit none

  integer :: stat_var = 0

  if (num_images() .lt. 3) error stop "Need at least three images to test."

  associate (me => this_image())
    if (me == 1) then
       continue !! image 1 does not participate and exits, creating a stopped image
    else
       associate (lhs => merge(me - 1, num_images(), me /= 1), &
              rhs => merge(me + 1, 1, me /= num_images()))
         sync images([lhs, rhs], STAT=stat_var)
         !! Only images bordering image 1 (i.e., 2 and `num_images()`) can
         !! accurately test whether a stopped image is present. All other
         !! images could be up ahead.
         if (stat_var /= STAT_STOPPED_IMAGE .and. me == 2) &
              error stop "Error: stat_var /= STAT_STOPPED_IMAGE: "
         if (stat_var /= STAT_STOPPED_IMAGE .and. me == num_images()) &
              error stop "Error: stat_var /= STAT_STOPPED_IMAGE: "
         if(me == 2) print *, 'Test passed.'
       end associate
    end if
  end associate
end program
