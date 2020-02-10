! This test checks if allocate on coarray variables acts as a barrier.

program alloc_as_barrier
  implicit none

  integer :: me

  call test_alloc(me)

  if(this_image() == 2) then
    if(me /= 1) then
      write(*,*) "Test failed.",me
    else
      write(*,*) "Test passed."
    endif
  endif

contains

  subroutine test_alloc(me)
    integer,intent(out) :: me
    integer,allocatable :: a(:)[:]

    me = this_image()
    if(me == 1) call sleep(1)
    allocate(a(10)[*],source=me)
    if(me > 1) me = a(2)[this_image()-1]
    deallocate(a)

  end subroutine

end program
