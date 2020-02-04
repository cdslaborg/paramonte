! SYNC IMAGES(*) with the STAT=STAT_STOPPED_IMAGE specifier
! Based on a test taken from UH caf-testsuite

program sync_images_stat
  use, intrinsic:: iso_fortran_env
  implicit none

  integer :: stat_var = 0, me

  me = this_image()

  if (me /= 1 ) then
     sync images(*,STAT=stat_var)
     if ( stat_var /= STAT_STOPPED_IMAGE) then
        print *, "Error:stat_var /= STAT_STOPPED_IMAGE: ", me
        ERROR STOP 1
     end if
     if(me == 2) print *, 'Test passed.'
  end if

  ! Image 1 implicitly synchronizes as part of normal termination
end program sync_images_stat
