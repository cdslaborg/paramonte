program main
   !! OpenCoarrays issue #488 MWE
   !! ===========================
   !!
   !! https://github.com/sourceryinstitute/OpenCoarrays/issues/488
   !!
   !! Running with 8 images, scattering from rank 1 fails when using a coarray
   !! in a derived type, whereas we observe no issue with pure coarrays.
   !!
   !! Should be run with 8 images as follows:
   !!    cafrun -np 8 ./a.out false
   !!    cafrun -np 8 ./a.out true
   implicit none

   type co_arr
      real, allocatable :: a(:, :)[:, :]
   end type

   type(co_arr) :: co1 ! derived type with coarray
   real, allocatable :: co2(:, :)[:, :] ! pure coarray
   real, allocatable :: glob(:, :), buf2d(:, :)
   integer, parameter :: nx = 2, ny = 3 ! local chunk, same over all images
   integer :: me, nimg, mx, my, i, j, k
   character(len = *), parameter :: tab = char(9), nl = char(10)
   character(len = 10) :: arg
   logical :: switch, test_passed

   test_passed = .false.
   call get_command_argument(1, arg)
   read (arg, '(L)') switch

   me = this_image()
   nimg = num_images()
   mx = nimg / 2 ! 2d grid distribution
   my = nimg / mx

   if (nimg /= 8) then
      stop 'example for 8 images'
   end if

   allocate(co1 % a(nx, ny)[mx, *])
   allocate(co2(nx, ny)[mx, *])
   allocate(buf2d(nx, ny), glob(mx * nx, my * ny))

   if (me == 1) print *, 'global size', [(mx * nx), (my * ny)], nl, 'local size', [nx, ny]

   ! call random_number(glob)
   glob = reshape([(i, i = 1, (mx * nx) * (my * ny))], shape(glob))

   sync all ! separate segments
   if (me == 1) then
      ! scatter glob from root (1st image)

      ! local to local
      co1 % a = glob(:nx, :ny)
      co2 = glob(:nx, :ny)

      ! loop over all other images
      do i = 1, mx
         do j = 1, my
            k = image_index(co2, [i, j])
            if (k /= 1) then
               buf2d = glob((i - 1) * nx + 1:i * nx, (j - 1) * ny + 1:j * ny) ! send buffer
               ! print *, 'filling up image #', k, '[', [i, j], ']', nl, tab, '==>', buf2d, nl
               if (switch) then
                  co1 % a(:nx, :ny)[i, j] = buf2d ! <= failure
               else
                  co2(:nx, :ny)[i, j] = buf2d ! ok
               end if
            end if
         end do
      end do
   end if
   sync all

   ! output local arrays after scattering
   ! print *, me, 'co1', co1 % a, nl, me, 'co2', co2

   ! collect results on rank 1
   call co_sum(co1 % a, result_image = 1)
   call co_sum(co2, result_image = 1)

   sync all
   test_passed = abs(sum(glob) - sum(merge(co1 % a, co2, switch))) < epsilon(0.)
   if (me == 1) then
      print *, 'all close ?', test_passed
      if(test_passed) then
         write(*,*) 'Test passed.'
      else
         write(*,*) 'Test failed!'
      endif
   end if
   deallocate(co1 % a, co2, glob)
end program
