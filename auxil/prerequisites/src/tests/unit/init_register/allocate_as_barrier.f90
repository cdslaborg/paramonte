! This test checks if allocate on coarray variables acts as a barrier.

program alloc_as_barrier
  implicit none

  integer :: me[*]
  integer,allocatable :: a(:)[:]

  me = this_image()

  if(me == 1) call sleep(1)

  allocate(a(10)[*])

  if(me > 1) then
    a = me[me-1]
    if(any(a /= me-1)) then
      write(*,*) "Test failed."
    else
      write(*,*) "Test passed."
    endif
  endif

end program
