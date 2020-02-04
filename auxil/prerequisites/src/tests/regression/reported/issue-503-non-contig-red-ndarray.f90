program main
   real, allocatable, dimension(:, :, :, :) :: arr, arr_copy
   integer :: i, me, nimg
   real :: redc, sumc, sumnc, rednc
   real :: eps = epsilon(0.)
   logical :: failc, failnc

   allocate(arr(9, 20, 3, 18))

   me = this_image()
   nimg = num_images()

   if (me == 1) then
      arr = reshape([(i, i=1, size(arr))], shape(arr))
   end if

   call co_broadcast(arr, source_image=1)
   arr_copy = arr

   print *, '==> CONTIGUOUS <=='
   sumc = sum(arr)
   redc = sumc
   call co_sum(redc)
   print *, 'sumc=', sumc
   print *, 'excpected: nimg * sumc=', nimg * sumc
   print *, 'got: redc=', redc

   failc = abs(redc - nimg * sumc) > eps

   print *, '==> NON CONTIGUOUS <=='
   sumnc = sum(arr(::3, ::2, :, ::5))

   call co_sum(arr_copy)
   call co_sum(arr(::3, ::2, :, ::5))

   rednc = sum(arr(::3, ::2, :, ::5))
   print *, 'sumnc=', sumnc
   print *, 'expected (nimg * sumnc)=', nimg * sumnc
   print *, 'expected=', sum(arr_copy(::3, ::2, :, ::5))
   print *, 'got: rednc=', rednc

   failnc = abs(rednc - nimg * sumnc) > eps

   sync all

   if (failc .or. failnc) then
      write(*, *) 'Test failed!'
      error stop 5
   else
      write(*, *) 'Test passed.'
   end if
   deallocate(arr)
end program
