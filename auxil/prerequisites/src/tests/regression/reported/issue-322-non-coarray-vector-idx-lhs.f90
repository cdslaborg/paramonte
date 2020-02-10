program mod_vecsub_01
  !! Issue #322 reported by @reinh-bader
  !! https://github.com/sourceryinstitute/OpenCoarrays/issues/322
  !! Related to other vector indexing issues
  implicit none
  integer, parameter :: ndim = 5, vdim = 2
  real :: vec(ndim), res(ndim)[*]
  integer :: idx(vdim)
  integer :: i, me
  logical :: ok[*]

  res = [ (real(i), i=1, ndim) ]
  vec = 0.0
  idx = [ ndim, 1 ]
  ok = .true.
  sync all
  me = this_image()
  vec(idx) = res(1:2)[1]
  if (vec(1) /= 2.0 .or. vec(5) /= 1.0) then
    critical
      ok[1] = .false.
      write(*, *) 'FAIL on image ',me,vec(idx)
    end critical
  end if
  if (me == 1) then
     if (ok) then
       write(*, *) 'Test passed.'
     else
        error stop 1
     end if
  end if
end program
