program slice
   type coarr
      real, allocatable :: a(:, :, :)[:, :, :]
   end type

   type(coarr) :: co
   integer :: nimg, me, z, nx, ny, nz, north, south, mex, mey, mez, coords(3)
   integer :: shape2d(2), shape3d(3)
   real, allocatable :: buf3d(:, :, :) ! a plane slice as a rank 3 array with a single transverse layer
   real, allocatable :: buf2d(:, :) ! a plane (2d) slice, normal in the y direction

   nx = 6
   ny = 4
   nz = 2

   me = this_image()
   nimg = num_images()

   if (nimg /= 8) stop

   allocate(co % a(nx, ny, nz)[1:2, 1:2, *])
   allocate(buf2d(nx, nz), buf3d(nx, 1, nz))

   ! this example should NOT reallocate buf2d nor buf3d
   ! compare shapes before and after syncing
   shape2d = shape(buf2d)
   shape3d = shape(buf3d)

   co % a = reshape([(z, z=1, nx * ny * nz)], shape(co % a))

   coords = this_image(co % a)
   mex = coords(1)
   mey = coords(2)
   mez = coords(3)

   north = mey + 1
   south = mey - 1

   sync all
   if (north <= 2) then
      z = image_index(co % a, [mex, north, mez])
      sync images(z)
      ! no reduction on rank
      buf3d = co % a(1:nx, 1:1, 1:nz)[mex, north, mez]
      co % a(1:nx, ny:ny, 1:nz) = buf3d

      ! reduction along dim 2
      buf2d = co % a(1:nx, 1, 1:nz)[mex, north, mez]
      co % a(1:nx, ny, 1:nz) = buf2d
   end if
   if (south >= 1) then
      z = image_index(co % a, [mex, south, mez])
      sync images(z)
      buf3d = co % a(1:nx, ny:ny, 1:nz)[mex, south, mez]
      co % a(1:nx, 1:1, 1:nz) = buf3d

      buf2d = co % a(1:nx, ny, 1:nz)[mex, south, mez]
      co % a(1:nx, 1, 1:nz) = buf2d
   end if
   sync all

   deallocate(co % a, buf2d, buf3d)
   
   if (any(shape2d /= shape(buf2d)) .or. any(shape3d /= shape(buf3d))) then
      write(*, *) 'Test failed!'
      error stop 5
   else
      write(*, *) 'Test passed.'
   end if

   ! Regression would cause error message:
   ! Fortran runtime error on image <...>: libcaf_mpi::caf_get_by_ref(): rank out of range.
end program
