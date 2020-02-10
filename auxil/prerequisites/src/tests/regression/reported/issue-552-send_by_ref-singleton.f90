program main
   implicit none

   type co_obj ! the test does not fail with pure coarrays, only when using coarrays in derived types
      real, allocatable, dimension(:, :, :) :: a[:], b, c, d, e
   end type

   type(co_obj) :: co
   integer :: me, remote, nimg, i, j, k, ni, nj, nk
   logical :: fail

   me = this_image()
   nimg = num_images()
   remote = merge(2, 1, me == 1)

   if (nimg /= 2) error stop 1

   ni = 2
   nj = 8
   nk = 4
   allocate(co % a(ni, nj, nk)[*])
   allocate(co % b(ni, nj, nk))
   allocate(co % c, co % d, co % e, mold=co % b)
   call random_number(co % b) ! the answer is a random array
   
   sync all
   co % a = co % b

   sync all
   co % c(1, :, :) = co % a(1, :, :)[remote] ! getter

   sync all
   co % e(1, :, :) = dble(111111111111111_8 * me) / 10**8
   co % a(1, :, :) = co % e(1, :, :)

   sync all
   ! singleton on the 1st dimension (delta == 1), triggers the bug

   ! FIXME: NOK, test FAILS
   ! co % a(1, :, :)[remote] = co % b(1, :, :) ! setter, (CAF_ARR_REF_SINGLE, CAF_ARR_REF_FULL, CAF_ARR_REF_FULL)

   ! OK, test pass, with patch
   co % a(1:1, :, :)[remote] = co % b(1:1, :, :) ! setter, (CAF_ARR_REF_RANGE, CAF_ARR_REF_FULL, CAF_ARR_REF_FULL)

   ! OK, test pass, without patch
   ! co % a(:, :, :)[remote] = co % b(:, :, :) ! setter, (CAF_ARR_REF_FULL,) * 3

   sync all
   co % d(1, :, :) = co % a(1, :, :)[remote] ! check

   sync all
   ! sequential flush (more readable)
   if (me == 2) sync images(1)
   write(*, *) ': : (set from remote on myself) : :', me
   do j = 1, nj; write(*, *) co % a(1, j, :); end do
   write(*, *) ': : (neighbor answer) : :', me
   do j = 1, nj; write(*, *) co % c(1, j, :); end do
   write(*, *) ': : (my answer) : :', me
   do j = 1, nj; write(*, *) co % b(1, j, :); end do
   write(*, *) ': : (what I see on remote) : :', me
   do j = 1, nj; write(*, *) co % d(1, j, :); end do
   write(*, *) ': : (initial local coarray) : :', me
   do j = 1, nj; write(*, *) co % e(1, j, :); end do
   write(*, *) ': : : :', me
   if (me == 1) sync images(2)
   
   sync all
   
   fail = any(abs(co % b(1, :, :) - co % d(1, :, :)) > epsilon(0.))
   
   if (fail) then
      write(*, *) 'Test failed!'
      error stop 5
   else
      write(*, *) 'Test passed.'
   end if

end program
