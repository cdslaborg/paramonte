!
! This program does a correctness check for
! ARRAY = SCALAR and ARRAY = ARRAY
!
program main
  implicit none
  integer, parameter :: n = 3
  integer, parameter :: m = 4

  ! Allocatable coarrays
  call one(-5, 1)
  call one(0, 0)
  call one(1, -5)
  call one(0, -11)

  ! Static coarrays
  call two()
  call three()
  write(*,*) 'Test passed'
contains
  subroutine one(lb1, lb2)
    integer, value :: lb1, lb2

    integer :: i_sgn1, i_sgn2, i, i_e, i_s, j, j_e, j_s
    integer, allocatable :: caf(:,:)[:]
    integer, allocatable :: a(:,:), b(:,:)

    allocate(caf(lb1:n+lb1-1, lb2:m+lb2-1)[*], &
         a(lb1:n+lb1-1, lb2:m+lb2-1), &
         b(lb1:n+lb1-1, lb2:m+lb2-1))

    b = reshape([(i*33, i = 1, size(b))], shape(b))

    ! Whole array: ARRAY = SCALAR
    a = -42
    caf = -42
    if(this_image() == num_images()) then
       caf = b
    endif
    sync all
    if (this_image() == 1) then
      a(:,:) = caf(lb1,lb2)[num_images()]
      print *, this_image(), '//', a, '//', b(lb1,lb2)
      print *, '>>>', any(a /= b(lb1,lb2))
      if (any (a /= b(lb1,lb2))) then
! FIXME: ABORTS UNLESS THERE IS SOME OTHER CODE
print *, 'HELLO!!!!!!!!!!!!!!!!!'
        error stop
      end if
    end if

    ! Whole array: ARRAY = ARRAY
    a = -42
    caf = -42
    if(this_image() == num_images()) then
       caf = b
    endif
    sync all
    if (this_image() == 1) then
      a(:,:) = caf(:,:)[num_images()]
      if (any (a /= b)) &
!FIXME
        print *, a
        print *, b
        print *, 'WRONG:', any (a /= b)
        error stop
      end if
    end if

    ! Array sections with different ranges and pos/neg strides
    do i_sgn1 = -1, 1, 2
      do i_sgn2 = -1, 1, 2
        do i=lb1, n+lb1-1
          do i_e=lb1, n+lb1-1
            do i_s=1, n
              do j=lb2, m+lb2-1
                do j_e=lb2, m+lb2-1
                  do j_s=1, m
                    ! ARRAY = SCALAR
                    a = -12
       		    caf = -42
    		    if(this_image() == num_images()) then
       		       caf = b
    		    endif
                    sync all
                    if (this_image() == 1) then
		      b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) &
			   = caf(lb1, lb2)[num_images()]
                    end if
                    sync all

                    ! ARRAY = ARRAY
                    a = -12
		    b = -32
                    if(this_image() == num_images()) then
               	       caf = a
              	     else
                       caf = -42
               	    endif
                    sync all
                    if (this_image() == 1) then
!		      b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) &
!			   = caf(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)[num_images()]
                    end if
                    sync all

                    if (this_image() == 1) then
 !                     if (any (a /= b)) then
 !                       print '(*(g0))', "bounds: ", lb1,":",n+lb1-1,", ", &
 !                            lb2,":",m+lb2-1
 !                       print '(*(g0))', "section: ", i,":",i_e,":",i_s*i_sgn1, &
 !                            ", ", j,":",j_e,":",j_s*i_sgn2
 !                       print *, i
 !                       print *, a
 !                       print *, caf
 !                       print *, a-caf
 !                       error stop
 !                     endif
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine one

  subroutine two()
    integer, parameter :: lb1 = -5, lb2 = 1

    integer :: i_sgn1, i_sgn2, i, i_e, i_s, j, j_e, j_s
    integer, save :: caf(lb1:n+lb1-1, lb2:m+lb2-1)[*]
    integer, save :: a(lb1:n+lb1-1, lb2:m+lb2-1)
    integer, save :: b(lb1:n+lb1-1, lb2:m+lb2-1)

    b = reshape([(i*33, i = 1, size(b))], shape(b))

    ! Whole array: ARRAY = SCALAR
    a = -12
    b = -32
    if(this_image() == num_images()) then
      caf = a
    else
      caf = -42
    endif
    sync all
    if (this_image() == 1) then
      b(:,:) = caf(lb1,lb2)[num_images()]
    end if
    sync all
    if (this_image() == 1) then
      if (any (a /= b)) &
           error stop
    end if

    ! Whole array: ARRAY = ARRAY
    a = -12
    b = -32
    if(this_image() == num_images()) then
      caf = a
    else
      caf = -42
    endif
    sync all
    if (this_image() == 1) then
      b(:,:) = caf(:,:)[num_images()]
    end if
    sync all
    if (this_image() == 1) then
      if (any (a /= b)) &
           error stop
    end if

    ! Array sections with different ranges and pos/neg strides
    do i_sgn1 = -1, 1, 2
      do i_sgn2 = -1, 1, 2
        do i=lb1, n+lb1-1
          do i_e=lb1, n+lb1-1
            do i_s=1, n
              do j=lb2, m+lb2-1
                do j_e=lb2, m+lb2-1
                  do j_s=1, m
                    ! ARRAY = SCALAR
                     a = -12
                     b = -32
                     if(this_image() == num_images()) then
                        caf = a
                     else
                        caf = -42
                     endif
                    sync all
                    if (this_image() == 1) then
                       b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) = caf(lb1,lb2)[num_images()]
                    end if
                    sync all

                    ! ARRAY = ARRAY
                    b = -32
                    a = -12
                    if(this_image() == num_images()) then
                       caf = a
                    else
                       caf = -42
                    endif
                    sync all
                    if (this_image() == 1) then
!                       b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) &
!                            =caf(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)[num_images()]
                    end if
                    sync all

                    if (this_image() == 1) then
 !                     if (any (a /= b)) then
 !                       print '(*(g0))', "bounds: ", lb1,":",n+lb1-1,", ", &
 !                            lb2,":",m+lb2-1
 !                       print '(*(g0))', "section: ", i,":",i_e,":",i_s*i_sgn1, &
 !                            ", ", j,":",j_e,":",j_s*i_sgn2
 !                       print *, i
 !                       print *, a
 !                       print *, caf
 !                       print *, a-caf
 !                       error stop
 !                     endif
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine two

  subroutine three()
    integer, parameter :: lb1 = 0, lb2 = 0

    integer :: i_sgn1, i_sgn2, i, i_e, i_s, j, j_e, j_s
    integer, save :: caf(lb1:n+lb1-1, lb2:m+lb2-1)[*]
    integer, save :: a(lb1:n+lb1-1, lb2:m+lb2-1)
    integer, save :: b(lb1:n+lb1-1, lb2:m+lb2-1)

    b = reshape([(i*33, i = 1, size(b))], shape(b))

    ! Whole array: ARRAY = SCALAR
    a = -12
    b = -32
    if(this_image() == num_images()) then
      caf = a
    else
      caf = -42
    endif
    sync all
    if (this_image() == 1) then
       b(:,:) = caf(lb1,lb2)[num_images()]
    end if
    sync all
    if (this_image() == 1) then
      if (any (a /= b)) &
           error stop
    end if

    ! Whole array: ARRAY = ARRAY
    a = -12
    b = -32
    if(this_image() == num_images()) then
      caf = a
    else
      caf = -42
    endif
    sync all
    if (this_image() == 1) then
       b(:,:) = caf(:,:)[num_images()]
    end if
    sync all
    if (this_image() == 1) then
      if (any (a /= b)) &
           error stop
    end if

    ! Array sections with different ranges and pos/neg strides
    do i_sgn1 = -1, 1, 2
      do i_sgn2 = -1, 1, 2
        do i=lb1, n+lb1-1
          do i_e=lb1, n+lb1-1
            do i_s=1, n
              do j=lb2, m+lb2-1
                do j_e=lb2, m+lb2-1
                  do j_s=1, m
                    ! ARRAY = SCALAR
                     a = -12
                     b = -32
                     if(this_image() == num_images()) then
                        caf = a
                     else
                        caf = -42
                     endif
                    sync all
                    if (this_image() == 1) then
                       b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) &
                            = caf(lb1,lb2)[num_images()]
                    end if
                    sync all

                    ! ARRAY = ARRAY
                     a = -12
                     b = -32
                     if(this_image() == num_images()) then
                        caf = a
                     else
                        caf = -42
                     endif
                     sync all
                    if (this_image() == 1) then
!                       b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) &
!                            = caf(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)[num_images()]
                    end if
                    sync all

                    if (this_image() == 1) then
!                      if (any (a /= b)) then
!                        print '(*(g0))', "bounds: ", lb1,":",n+lb1-1,", ", &
!                             lb2,":",m+lb2-1
!                        print '(*(g0))', "section: ", i,":",i_e,":",i_s*i_sgn1, &
!                             ", ", j,":",j_e,":",j_s*i_sgn2
!                        print *, i
!                        print *, a
!                        print *, caf
!                        print *, a-caf
!                        error stop
!                      endif
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine three
end program main
