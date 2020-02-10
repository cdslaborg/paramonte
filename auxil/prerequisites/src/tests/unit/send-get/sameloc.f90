! This program tests the capability of copying data on the
! same memory location and within the same image from
! different memory locations.
program sameloc
  implicit none

  integer,codimension[*] :: a
  integer,dimension(10),codimension[*] :: b,c
  integer,dimension(9,10),codimension[*] :: m
  integer,dimension(10) :: t
  integer :: i,j
  logical :: tests_passed

  tests_passed = .true.

  a = 10
  b(1:5) = 1
  b(6:10) = -1
  c(1:5) = 1
  c(6:10) = -1

  t(:) = b(:)
  t(1:5) = b(2:6)

  do i=1,9
    m(i,:) = (/ (j, j = 1, 10) /)
  enddo

  sync all

  a = a[1]
  if (this_image() == 1) write(*,*) 'OK',a

  t = (/ (j, j = 1, 10) /)

  if(this_image() == 1) then
    c = m(1,:)[1]
    if(any(c(:) /= t(:))) then
      tests_passed = .false.
      error stop "get row failed"
    else
     tests_passed = tests_passed .and. .true.
     write(*,*) 'ok get row'
    endif
  endif

  sync all

  if(this_image() == 1) then
    do i=1,10
      if(m(9,i)[1] /= t(i)) then
	write(*,*) 'pos',i,'value get',m(9,i)[1],'value t',t(i)
        tests_passed = .false.
        error stop "get element from matrix failed"
      else
        tests_passed = tests_passed .and. .true.
      endif
    enddo
  endif

  if (this_image() == 1) write(*,*) 'Ok get element from matrix'

  sync all

  m(9,:) = 1

  if(this_image() == 1) then
    do i=1,10
      m(9,i)[1] = i
      if(m(9,i)[1] /= t(i)) then
        write(*,*) 'pos',i,'value get',m(9,i)[1],'value t',t(i)
        tests_passed = .false.
        error stop "put element from matrix failed"
      else
        tests_passed = tests_passed .and. .true.
      endif
    enddo
  endif

  if (this_image() == 1) write(*,*) 'Ok put element from matrix'

  t(:) = b(:)
  t(1:5) = b(2:6)

  c(1:5) = 1
  c(6:10) = -1

  sync all

  if(this_image() == 1) then
    b(1:5)[1] = b(2:6)
    if(any(b(:) /= t(:))) then
      tests_passed = .false.
      error stop "put overlapped failed"
    else
      tests_passed = tests_passed .and. .true.
      write(*,*) 'OK put overlapped'
    endif
  endif

  b(1:5) = 1
  b(6:10) = -1

  sync all

  if(this_image() == 1) then
    b(1:5)[1] = b(2:6)[1]
    if(any(b(:) /= t(:))) then
      tests_passed = .false.
      error stop "putget overlapped failed"
    else
      tests_passed = tests_passed .and. .true.
      write(*,*) 'OK putget overlapped'
    endif
  endif

  t(:) = c(:)
  t(10:1:-1) = t(:)

  sync all

  if(this_image() == 1) then
    c(10:1:-1)[1] = c(:)
    if(any(t(:) /= c(:))) then
      tests_passed = .false.
      write(*,*) 'Error in put reversed'
      write(*,*) c
      write(*,*) t
      error stop "put reversed failed"
    else
      tests_passed = tests_passed .and. .true.
      write(*,*) 'OK put reversed'
    endif
  endif

  c(1:5) = 1
  c(6:10) = -1

  t(:) = c(:)
  t(10:1:-1) = t(:)

  if(this_image() == 1) then
    c(:) = c(10:1:-1)[1]
    if(any(t(:) /= c(:))) then
      tests_passed = .false.
      write(*,*) c
      write(*,*) t
      error stop "get reversed failed"
    else
      tests_passed = tests_passed .and. .true.
      write(*,*) 'OK get reversed'
    endif
  endif

  if ( .not. tests_passed ) then
     error stop "Test failures exist!"
  end if

  sync all

  if ( tests_passed ) then
    if (this_image() == 1) write(*,*) 'Test passed'
  end if

end program
