program source_allocation_no_sync
  !! author: Damian Rouson and Izaak Beekman
  !! category: regression
  !!
  !! [Issue #243](https://github.com/sourceryinstitute/opencoarrays/issues/243)
  !!
  !! [GFortran PR 78505](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=78505):
  !!
  !! @note The test must be run with less than or equal to 32 images,
  !! and the number of images must be a power of two. Valid numbers of
  !! images are: 2, 4, 8, 16, or 32
  !!
  !! Sourced allocation of a coarray object performs a synchronization
  !! after the allocation, BUT *before* the assignment of
  !! `source`. This violates the
  !! [standard](http://open-std.org/JTC1/SC22/WG5/5559) Section
  !! 9.7.1.2, paragraph 4:
  !!
  !! > When an ALLOCATE statement is executed for which an
  !! > allocate-object is a coarray, there is an implicit
  !! > synchronization of all active images in the current team. On
  !! > those images, if no error condition other than
  !! > STAT_STOPPED_IMAGE or STAT_FAILED_IMAGE occurs, execution of
  !! > the segment (11.6.2) following the statement is delayed until
  !! > all other active images in the current team have executed the
  !! > same statement the same number of times in this team. The
  !! > coarray shall not become allocated on an image unless it is
  !! > successfully allocated on all active images in this team.
  !!

  implicit none
  integer, allocatable :: f(:)[:]
  integer, parameter :: num_points=32
  integer :: me,ni,my_num_points,neighbor_last_element
  me = this_image()
  if (mod(num_points,num_images())/=0) error stop "num_points not evenly divisible by num_images()"
  my_num_points = num_points/num_images()
  allocate( f(my_num_points)[*], source = 1 )
  if (me>1) then
     neighbor_last_element = f(my_num_points)[me-1]
     if (neighbor_last_element /=1) then
        print *,"Image ",me," gets ",neighbor_last_element
        error stop "Synchronization did not happen after assignment in sourced allocation!"
     end if 
  end if
  sync all
  if ( me == 1 ) then
     write(*,'(a)') "Test passed."
  end if
end program source_allocation_no_sync
