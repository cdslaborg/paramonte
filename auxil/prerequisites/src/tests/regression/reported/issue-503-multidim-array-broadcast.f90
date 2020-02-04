! caf -o bcast bcast.f90
! cafrun -np 2 ./bcast
! mwe for issue github.com/sourceryinstitute/OpenCoarrays/issues/503

program main
   real, allocatable, dimension(:, :, :) :: arr1, arr2
   integer :: i, me, nimg
   real :: red1, sum1, red2, sum2
   
   allocate(arr1(10, 20, 8))
   allocate(arr2(10, 20, 30))

   me = this_image()
   nimg = num_images()

   if (me == 1) then
      arr1 = reshape([(i, i=1, size(arr1))], shape(arr1))
      arr2 = reshape([(i, i=1, size(arr2))], shape(arr2))
   end if

   call co_broadcast(arr1, source_image=1)
   call co_broadcast(arr2, source_image=1)

   sum1 = sum(arr1)
   sum2 = sum(arr2)

   red1 = sum1
   red2 = sum2

   call co_sum(red1)
   call co_sum(red2)

   sync all
   print *, me, ' sum1=', sum1, ' red1=', red1
   print *, me, ' sum2=', sum2, ' red2=', red2

   if (abs(red1 - nimg * sum1) > epsilon(0.) .or. abs(red2 - nimg * sum2) > epsilon(0.)) then
      write(*,*) 'Test failed!'
      error stop 5
   else
      write(*,*) 'Test passed.'
   end if
   deallocate(arr1, arr2)
end program
