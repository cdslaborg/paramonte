! caf -o gather gather.f90
! cafrun -np 2 ./gather
! this example mimics a mpi_gatherv with root=1, and variable size chunks
! gather into specified locations from all processes in a group

program gather
   implicit none

   type gvec ! a global vector
      real, allocatable :: a(:)
   end type

   type(gvec) :: gc[*]
   real, allocatable :: gv(:), tmp(:)
   integer :: me, nimg, gsize, i, lo, hi
   logical :: fail = .false.

   me = this_image()
   nimg = num_images()

   allocate(gc % a(2 * me)) ! variable size data in container
   gc % a = [(me * i, i=1, 2 * me)] ! assignement

   ! collect the global vector size by summing local sizes
   gsize = size(gc % a)
   call co_sum(gsize, result_image=1)
   sync all

   if (me == 1) then
      if (gsize /= 6) error stop 1
      allocate(gv(gsize)) ! allocate a global vector of size 6 on img 1
      lo = 1
      do i = 1, nimg
         tmp = gc[i] % a ! note: automatic reallocation of tmp
         hi = lo + size(tmp) - 1
         gv(lo:hi) = tmp
         lo = hi + 1 ! start of next chunk
      end do
      print *, 'gv=', gv, ' sum=', sum(gv)

      if (abs(sum(gv) - 23.) > epsilon(0.)) fail = .true.
   end if

   sync all

   ! CMake test output handler
   call co_broadcast(fail, source_image=1)
   if (fail) then
      write(*, *) 'Test failed!'
      error stop 5
   else
      write(*, *) 'Test passed.'
   end if

end program
