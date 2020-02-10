program atomic
use iso_fortran_env
implicit none

integer :: me,np,res
integer(atomic_int_kind) :: atom[*]

me = this_image()
np = num_images()

call atomic_define(atom[1],0)

sync all

call ATOMIC_ADD (atom[1], me)

sync all

if(me == 1) then
  call atomic_ref(res,atom[1])
  if(res /= (np*(np+1))/2) then
    write(*,*) 'res',res
    error stop "Atomic ref failed"
  endif
  write(*,*) 'OK'
endif

sync all

if (me == 1) then
   write(*,*) "Test passed"
end if

end program
