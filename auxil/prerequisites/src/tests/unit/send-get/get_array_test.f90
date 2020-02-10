!
! This program does a correctness check for
! ... = ARRAY[idx] and ... = SCALAR[idx]
!
program main
  implicit none
  integer, parameter :: n = 3
  integer, parameter :: m = 4

  logical :: error_printed=.false.

  ! Allocatable coarrays
  call one(-5, 1)
  call one(0, 0)
  call one(1, -5)
  call one(0, -11)

  ! Static coarrays
  call two()
  call three()

  if (error_printed)  error stop
  sync all

  if (this_image()==1) print *,'Test passed.'
contains
  subroutine one(lb1, lb2)
    integer, value :: lb1, lb2

    integer :: i_sgn1, i_sgn2, i, i_e, i_s, j, j_e, j_s
    integer, allocatable :: caf(:,:)[:]
    integer, allocatable :: a(:,:), b(:,:), c(:,:)

    allocate(caf(lb1:n+lb1-1, lb2:m+lb2-1)[*], &
         a(lb1:n+lb1-1, lb2:m+lb2-1), &
         b(lb1:n+lb1-1, lb2:m+lb2-1), &
         c(lb1:n+lb1-1, lb2:m+lb2-1))

    b = reshape([(i*33, i = 1, size(b))], shape(b))

    ! Whole array: ARRAY = SCALAR
    a = b
    caf = -42
    c = caf
    sync all
    if(this_image() == 1) then
       a(:,:) = caf(lb1,lb2)[num_images()]
    end if
    sync all
    if(this_image()==1) then
       if(any (a /= c)) call print_and_register( "ARRAY = SCALAR failed in get_array_test")
    endif

    ! Whole array: ARRAY = ARRAY
    caf = -42
    a = b
    c = caf
    if (this_image() == 1) then
       a(:,:) = caf(:,:)[num_images()]
    endif
    sync all
    if(this_image()==1) then
       if (any (a /= c)) then
          print *, 'RES 1:', any (a /= c)
          print *, a
          print *, c
          ! FIXME: Without the print lines above, it always fails. Why?
          call print_and_register( "ARRAY = ARRAY failed in get_array_test")
       end if
    endif

    ! Scalar assignment
    a = -42
    caf = -42
    c = caf
    sync all
    do j = lb2, m+lb2-1
       do i = n+lb1-1, lb1, -2
          a(i,j) = b(i,j)
       end do
    end do
    do j = lb2, m+lb2-1
       do i = lb1, n+lb1-1, 2
          a(i,j) = b(i,j)
       end do
    end do
    sync all
    if(this_image() == 1) then
       do j = lb2, m+lb2-1
          do i = n+lb1-1, lb1, -2
             a(i,j) = caf(i,j)[num_images()]
          end do
       end do
       do j = lb2, m+lb2-1
          do i = lb1, n+lb1-1, 2
             a(i,j) = caf(i,j)[num_images()]
          end do
       end do
    endif
    sync all
    if(this_image() == 1) then
       if (any (a /= c)) then
          print *, 'RES 2:', any (a /= c)
          print *, this_image(), ': ', a
          print *, this_image(), ': ', c
          ! FIXME: Without the print lines above, it always fails. Why?
          call print_and_register( "scalar assignment failed in get_array_test")
       end if
    endif
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
                            a = -42
                            caf = -42
                            c = a
                            a(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) = b(lb1, lb2)
                            sync all
                            if (this_image() == 1) then
                               a(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) = caf(lb1,lb2)[num_images()]
                            end if
                            sync all
                            if (this_image() == 1) then
                               if (any (a /= c)) then
                                  print '(*(g0))', "bounds: ", lb1,":",n+lb1-1,", ", &
                                       lb2,":",m+lb2-1
                                  print '(*(g0))', "section: ", i,":",i_e,":",i_s*i_sgn1, &
                                       ", ", j,":",j_e,":",j_s*i_sgn2
                                  print *, i
                                  print *, a
                                  print *, c
                                  print *, a-c
                                  call print_and_register( "array sections with ranges and strides failed in get_array_test")
                               endif
                            end if
                            ! ARRAY = ARRAY
                            caf = -42
                            a = -42
                            c = a
                            a(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) &
                                 = b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)
                            sync all
                            if (this_image() == 1) then
                               a(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) = caf(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)[num_images()]
                            end if
                            sync all

                            if (this_image() == 1) then
                               if (any (a /= c)) then
                                  print '(*(g0))', "bounds: ", lb1,":",n+lb1-1,", ", &
                                       lb2,":",m+lb2-1
                                  print '(*(g0))', "section: ", i,":",i_e,":",i_s*i_sgn1, &
                                       ", ", j,":",j_e,":",j_s*i_sgn2
                                  print *, i
                                  print *, a
                                  print *, c
                                  print *, a-c
                                  call print_and_register( "array sections with ranges and strides failed in get_array_test")
                               endif
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
    caf = -42
    a = -42
    a(:,:) = b(lb1, lb2)
    sync all
    if (this_image() == 1) then
      caf(:,:)[num_images()] = b(lb1, lb2)
    end if
    sync all
    if (this_image() == num_images()) then
      if (any (a /= caf)) &
           call print_and_register( "Array = scalar failed in subroutine two get_array_test")
    end if

    ! Whole array: ARRAY = ARRAY
    caf = -42
    a = -42
    a(:,:) = b(:, :)
    sync all
    if (this_image() == 1) then
      caf(:,:)[num_images()] = b(:, :)
    end if
    sync all
    if (this_image() == num_images()) then
      if (any (a /= caf)) &
           call print_and_register( "Array = array failed in subroutine two get_array_test")
    end if

    ! Scalar assignment
    caf = -42
    a = -42
    do j = lb2, m+lb2-1
      do i = n+lb1-1, 1, -2
        a(i,j) = b(i,j)
      end do
    end do
    do j = lb2, m+lb2-1
      do i = 1, n+lb1-1, 2
        a(i,j) = b(i,j)
      end do
    end do
    sync all
    if (this_image() == 1) then
      do j = lb2, m+lb2-1
        do i = n+lb1-1, 1, -2
          caf(i,j)[num_images()] = b(i, j)
        end do
      end do
      do j = lb2, m+lb2-1
        do i = 1, n+lb1-1, 2
          caf(i,j)[num_images()] = b(i, j)
        end do
      end do
    end if
    sync all
    if (this_image() == num_images()) then
      if (any (a /= caf)) &
           call print_and_register( "scalar assignment failed in subroutine two get_array_test")
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
                    caf = -42
                    a = -42
                    a(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) = b(lb1, lb2)
                    sync all
                    if (this_image() == 1) then
                      caf(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)[num_images()] &
                           = b(lb1, lb2)
                    end if
                    sync all

                    ! ARRAY = ARRAY
                    caf = -42
                    a = -42
                    a(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) &
                         = b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)
                    sync all
                    if (this_image() == 1) then
                      caf(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)[num_images()] &
                           = b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)
                    end if
                    sync all

                    if (this_image() == num_images()) then
                      if (any (a /= caf)) then
                        print '(*(g0))', "bounds: ", lb1,":",n+lb1-1,", ", &
                             lb2,":",m+lb2-1
                        print '(*(g0))', "section: ", i,":",i_e,":",i_s*i_sgn1, &
                             ", ", j,":",j_e,":",j_s*i_sgn2
                        print *, i
                        print *, a
                        print *, caf
                        print *, a-caf
                        call print_and_register( "arrays with ranges and strides failed sub. two get_array_test failed")
                      endif
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
    caf = -42
    a = -42
    a(:,:) = b(lb1, lb2)
    sync all
    if (this_image() == 1) then
      caf(:,:)[num_images()] = b(lb1, lb2)
    end if
    sync all
    if (this_image() == num_images()) then
      if (any (a /= caf)) &
           call print_and_register( "Array = scalar subroutine three get_array_test failed")
    end if

    ! Whole array: ARRAY = ARRAY
    caf = -42
    a = -42
    a(:,:) = b(:, :)
    sync all
    if (this_image() == 1) then
      caf(:,:)[num_images()] = b(:, :)
    end if
    sync all
    if (this_image() == num_images()) then
      if (any (a /= caf)) &
           call print_and_register( "Array = array subroutine three get_array_test failed")
    end if

    ! Scalar assignment
    caf = -42
    a = -42
    do j = lb2, m+lb2-1
      do i = n+lb1-1, 1, -2
        a(i,j) = b(i,j)
      end do
    end do
    do j = lb2, m+lb2-1
      do i = 1, n+lb1-1, 2
        a(i,j) = b(i,j)
      end do
    end do
    sync all
    if (this_image() == 1) then
      do j = lb2, m+lb2-1
        do i = n+lb1-1, 1, -2
          caf(i,j)[num_images()] = b(i, j)
        end do
      end do
      do j = lb2, m+lb2-1
        do i = 1, n+lb1-1, 2
          caf(i,j)[num_images()] = b(i, j)
        end do
      end do
    end if
    sync all
    if (this_image() == num_images()) then
      if (any (a /= caf)) &
           call print_and_register( "scalar assignment subroutine three get_array_test failed")
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
                    caf = -42
                    a = -42
                    a(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) = b(lb1, lb2)
                    sync all
                    if (this_image() == 1) then
                      caf(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)[num_images()] &
                           = b(lb1, lb2)
                    end if
                    sync all

                    ! ARRAY = ARRAY
                    caf = -42
                    a = -42
                    a(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2) &
                         = b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)
                    sync all
                    if (this_image() == 1) then
                      caf(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)[num_images()] &
                           = b(i:i_e:i_s*i_sgn1, j:j_e:j_s*i_sgn2)
                    end if
                    sync all

                    if (this_image() == num_images()) then
                      if (any (a /= caf)) then
                        print '(*(g0))', "bounds: ", lb1,":",n+lb1-1,", ", &
                             lb2,":",m+lb2-1
                        print '(*(g0))', "section: ", i,":",i_e,":",i_s*i_sgn1, &
                             ", ", j,":",j_e,":",j_s*i_sgn2
                        print *, i
                        print *, a
                        print *, caf
                        print *, a-caf
                        call print_and_register( "range stride in subroutine three get_array_test failed")
                      endif
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

  subroutine print_and_register(error_message)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_message
    write(error_unit,*) error_message
    error_printed=.true.
  end subroutine

end program main
