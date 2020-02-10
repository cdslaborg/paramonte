!****p* doc/sample-fortran-source.f90
! NAME
!   sendrcv
! SYNOPSIS
!   In this simple coarray Fortran program, image 1 puts its
!   local elements of an array coarray into the corresponding
!   elements of image 2.  The corresponding C program that an
!   OpenCoarrays-compatible compiler might generate from this
!   code is in doc/sample-compiler-output.c.
! SOURCE
program sendrecv
  use iso_fortran_env
  implicit none

  integer :: me, np, i
  integer, parameter :: n=1000
  real(kind=real64), allocatable :: d(:)[:]

  allocate(d(n)[*])

  np = num_images()
  me = this_image()

  do i=1,n
     d(i) = i
  enddo

  sync all

  if(me == 1) d(:)[2] = d

  sync all

  deallocate(d)

end program
!******
