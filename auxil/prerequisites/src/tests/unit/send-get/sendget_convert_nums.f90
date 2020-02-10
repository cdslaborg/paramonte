!! Thoroughly test sendget, i.e. foo [M] = bar[N] in all variants
!!
!! Do simple tests for sendget(). These test comprise
!!
!! FOO [M] = BAR [N]
!!
!! where
!!
!!  FOO                BAR                images
!! scalar            scalar              N == M == me
!!  int(k e [1,4])    int(k e [1,4])
!!  real(k e [4,8])   real(k e [4,8])
!!  int(k e [1,4])    real(k e [4,8])
!!  real(k e [4,8])   int(k e [1,4])
!!
!! array(1:5)        scalar
!!  int(k e [1,4])    int(k e [1,4])
!!  real(k e [4,8])   real(k e [4,8])
!!  int(k e [1,4])    real(k e [4,8])
!!  real(k e [4,8])   int(k e [1,4])
!!
!! array(1:5)        array(1:5)
!!  int(k e [1,4])    int(k e [1,4])
!!  real(k e [4,8])   real(k e [4,8])
!!  int(k e [1,4])    real(k e [4,8])
!!  real(k e [4,8])   int(k e [1,4])
!!
!! array(1:3)       array(::2)
!!  int(k e [1,4])    int(k e [1,4])
!!  real(k e [4,8])   real(k e [4,8])
!!  int(k e [1,4])    real(k e [4,8])
!!  real(k e [4,8])   int(k e [1,4])
!!
!! array(4:5)       array(2::2)
!!  int(k e [1,4])    int(k e [1,4])
!!  real(k e [4,8])   real(k e [4,8])
!!  int(k e [1,4])    real(k e [4,8])
!!  real(k e [4,8])   int(k e [1,4])
!!
!! array(1:3)      array(3:1:-1)
!!  int(k e [1,4])    int(k e [1,4])
!!  real(k e [4,8])   real(k e [4,8])
!!  int(k e [1,4])    real(k e [4,8])
!!  real(k e [4,8])   int(k e [1,4])
!!
!! all of the above but for             N != M != me
!!
!! And may be some other, I've forgotten.
!!
!! Author: Andre Vehreschild, 2017

program sendget_convert_nums

  implicit none

  real(kind=4), parameter :: tolerance4 = 1.0e-4
  real(kind=8), parameter :: tolerance4to8 = 1.0E-4
  real(kind=8), parameter :: tolerance8 = 1.0E-6

  integer(kind=1), codimension[*] :: co_int_scal_src_k1, co_int_scal_dst_k1
  integer(kind=4), codimension[*] :: co_int_scal_src_k4, co_int_scal_dst_k4
  real(kind=4)   , codimension[*] :: co_real_scal_src_k4, co_real_scal_dst_k4
  real(kind=8)   , codimension[*] :: co_real_scal_src_k8, co_real_scal_dst_k8

  integer(kind=1), dimension(1:5), codimension[*] :: co_int_src_k1, co_int_dst_k1
  integer(kind=4), dimension(1:5), codimension[*] :: co_int_src_k4, co_int_dst_k4
  real(kind=4)   , dimension(1:5), codimension[*] :: co_real_src_k4, co_real_dst_k4
  real(kind=8)   , dimension(1:5), codimension[*] :: co_real_src_k8, co_real_dst_k8

  logical :: error_printed=.false.

  associate(me => this_image(), np => num_images())
    if (np < 3) error stop 'Can not run with less than 3 images.'

    co_int_scal_src_k1 = 42
    co_int_scal_dst_k1 = -1
    co_int_scal_src_k4 = 42
    co_int_scal_dst_k4 = -1
    co_int_src_k1 = INT([ 5, 4, 3, 2, 1], 1)
    co_int_dst_k1 = -1
    co_int_src_k4 = [ 5, 4, 3, 2, 1]
    co_int_dst_k4 = -1

    co_real_scal_src_k4 = 37.042
    co_real_scal_src_k8 = REAL(37.042, 8)
    co_real_scal_dst_k4 = -1.0
    co_real_scal_dst_k8 = -1.0
    co_real_src_k4 = [ 5.1, 4.2, 3.3, 2.4, 1.5]
    co_real_src_k8 = REAL([ 5.1, 4.2, 3.3, 2.4, 1.5], 8)
    co_real_dst_k4 = -1.0
    co_real_dst_k8 = -1.0
    ! First check send/copy to self
    if (me == 1) then
      co_int_scal_dst_k1[1] = co_int_scal_src_k1[1]
      print *, co_int_scal_dst_k1
      if (co_int_scal_dst_k1 /= co_int_scal_src_k1) &
        call print_and_register( 'sendgetget scalar int kind=1 to kind=1 self failed.')

      co_int_scal_dst_k4[1] = co_int_scal_src_k4[1]
      print *, co_int_scal_dst_k4
      if (co_int_scal_dst_k4 /= co_int_scal_src_k4) &
        call print_and_register( 'sendgetget scalar int kind=4 to kind=4 self failed.')

      co_int_scal_dst_k4[1] = co_int_scal_src_k1[1]
      print *, co_int_scal_dst_k4
      if (co_int_scal_dst_k4 /= co_int_scal_src_k4) call print_and_register( 'sendget scalar int kind=1 to kind=4 self failed.')

      co_int_scal_dst_k1[1] = co_int_scal_src_k4[1]
      print *, co_int_scal_src_k1
      if (co_int_scal_dst_k1 /= co_int_scal_src_k1) call print_and_register( 'sendget scalar int kind=4 to kind=1 self failed.')

      co_int_dst_k1(:)[1] = co_int_src_k1(:)[1]
      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= co_int_src_k1)) call print_and_register( 'sendget int kind=1 to kind=1 self failed.')

      co_int_dst_k4(:)[1] = co_int_src_k4(:)[1]
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= co_int_src_k4)) call print_and_register( 'sendget int kind=4 to kind=4 self failed.')

      co_int_dst_k4(:)[1] = co_int_src_k1(:)[1]
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= co_int_src_k4)) call print_and_register( 'sendget int kind=1 to kind=4 self failed.')

      co_int_dst_k1(:)[1] = co_int_src_k4(:)[1]
      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= co_int_src_k1)) call print_and_register( 'sendget int kind=4 to kind=1 self failed.')
    else if (me == 2) then ! Do the real copy to self checks on image 2
      co_real_scal_dst_k4[2] = co_real_scal_src_k4[2]
      print *, co_real_scal_dst_k4
      if (abs(co_real_scal_dst_k4 - co_real_scal_src_k4) > tolerance4) &
        & call print_and_register( 'sendget scalar real kind=4 to kind=4 self failed.')

      co_real_scal_dst_k8[2] = co_real_scal_src_k8[2]
      print *, co_real_scal_dst_k8
      if (abs(co_real_scal_dst_k8 - co_real_scal_src_k8) > tolerance8) &
        & call print_and_register( 'sendget scalar real kind=8 to kind=8 self failed.')

      co_real_scal_dst_k8[2] = co_real_scal_src_k4[2]
      print *, co_real_scal_dst_k8
      if (abs(co_real_scal_dst_k8 - co_real_scal_src_k8) > tolerance4to8) &
        & call print_and_register( 'sendget scalar real kind=4 to kind=8 self failed.')

      co_real_scal_dst_k4[2] = co_real_scal_src_k8[2]
      print *, co_real_scal_dst_k4
      if (abs(co_real_scal_dst_k4 - co_real_scal_src_k4) > tolerance4) &
        & call print_and_register( 'sendget scalar real kind=8 to kind=4 self failed.')

      co_real_dst_k4(:)[2] = co_real_src_k4(:)[2]
      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - co_real_src_k4) > tolerance4)) &
        call print_and_register( 'sendget real kind=4 to kind=4 self failed.')

      co_real_dst_k8(:)[2] = co_real_src_k8(:)[2]
      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - co_real_src_k8) > tolerance8)) &
        call print_and_register( 'sendget real kind=8 to kind=8 self failed.')

      co_real_dst_k8(:)[2] = co_real_src_k4(:)[2]
      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - co_real_src_k8) > tolerance4to8)) &
        call print_and_register( 'sendget real kind=4 to kind=8 self failed.')

      co_real_dst_k4(:)[2] = co_real_src_k8(:)[2]
      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - co_real_src_k4) > tolerance4)) &
        call print_and_register( 'sendget real kind=8 to kind=4 self failed.')
    end if

    sync all
    if (me == 2) then
      co_int_scal_dst_k1[1] = co_int_scal_src_k1[3]

      co_int_scal_dst_k4[1] = co_int_scal_src_k4[3]

      co_int_dst_k1(:)[1] = co_int_src_k1(:)[3]

      co_int_dst_k4(:)[1] = co_int_src_k4(:)[3]

      co_real_scal_dst_k4[1] = co_real_scal_src_k4[3]

      co_real_scal_dst_k8[1] = co_real_scal_src_k8[3]

      co_real_dst_k4(:)[1] = co_real_src_k4(:)[3]

      co_real_dst_k8(:)[1] = co_real_src_k8(:)[3]
    end if

    sync all
    if (me == 1) then
      print *, co_int_scal_dst_k1
      if (co_int_scal_dst_k1 /= co_int_scal_src_k1) &
        call print_and_register( 'sendget scalar int kind=1 to kind=1 to image 2 failed.')

      print *, co_int_scal_dst_k4
      if (co_int_scal_dst_k4 /= co_int_scal_src_k4) &
        call print_and_register( 'sendget scalar int kind=4 to kind=4 to image 2 failed.')

      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= co_int_src_k1)) call print_and_register( 'sendget int kind=1 to kind=1 to image 2 failed.')

      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= co_int_src_k4)) call print_and_register( 'sendget int kind=4 to kind=4 to image 2 failed.')

      print *, co_real_scal_dst_k4
      if (abs(co_real_scal_dst_k4 - co_real_scal_src_k4) > tolerance4) &
        & call print_and_register( 'sendget scalar real kind=4 to kind=4 to image 2 failed.')

      print *, co_real_scal_dst_k8
      if (abs(co_real_scal_dst_k8 - co_real_scal_src_k8) > tolerance8) &
        & call print_and_register( 'sendget scalar real kind=8 to kind=8 to image 2 failed.')

      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - co_real_src_k4) > tolerance4)) &
        call print_and_register( 'sendget real kind=4 to kind=4 to image 2 failed.')

      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - co_real_src_k8) > tolerance8)) &
        call print_and_register( 'sendget real kind=8 to kind=8 to image 2 failed.')
    end if

    sync all
    if (me == 2) then
      co_int_scal_dst_k4[1] = co_int_scal_src_k1[3]

      co_int_scal_dst_k1[1] = co_int_scal_src_k4[3]

      co_int_dst_k4(:)[1] = co_int_src_k1(:)[3]

      co_int_dst_k1(:)[1] = co_int_src_k4(:)[3]

      co_real_scal_dst_k8[1] = co_real_scal_src_k4[3]

      co_real_scal_dst_k4[1] = co_real_scal_src_k8[3]

      co_real_dst_k8(:)[1] = co_real_src_k4(:)[3]

      co_real_dst_k4(:)[1] = co_real_src_k8(:)[3]
    end if

    sync all
    if (me == 1) then
      print *, co_int_scal_dst_k4
      if (co_int_scal_dst_k4 /= co_int_scal_src_k4) &
        call print_and_register( 'sendget scalar int kind=1 to kind=4 to image 2 failed.')

      print *, co_int_scal_dst_k1
      if (co_int_scal_dst_k1 /= co_int_scal_src_k1) &
        call print_and_register( 'sendget scalar int kind=4 to kind=1 to image 2 failed.')

      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= co_int_src_k4)) call print_and_register( 'sendget int kind=1 to kind=4 to image 2 failed.')

      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= co_int_src_k1)) call print_and_register( 'sendget int kind=4 to kind=1 to image 2 failed.')

      print *, co_real_scal_dst_k8
      if (abs(co_real_scal_dst_k8 - co_real_scal_src_k8) > tolerance4to8) &
        & call print_and_register( 'sendget scalar real kind=4 to kind=8 to image 2 failed.')

      print *, co_real_scal_dst_k4
      if (abs(co_real_scal_dst_k4 - co_real_scal_src_k4) > tolerance4) &
        & call print_and_register( 'sendget scalar real kind=8 to kind=4 to image 2 failed.')

      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - co_real_src_k8) > tolerance4to8)) &
        call print_and_register( 'sendget real kind=4 to kind=8 to image 2 failed.')

      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - co_real_src_k4) > tolerance4)) &
        call print_and_register( 'sendget real kind=8 to kind=4 to image 2 failed.')
    end if

    ! Scalar to array replication
    sync all
    if (me == 2) then
      co_int_dst_k4(:)[1] = co_int_scal_src_k4[3]

      co_int_dst_k1(:)[1] = co_int_scal_src_k1[3]

      co_real_dst_k8(:)[1] = co_real_scal_src_k8[3]

      co_real_dst_k4(:)[1] = co_real_scal_src_k4[3]
    end if

    sync all
    if (me == 1) then
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= co_int_scal_src_k4)) &
        call print_and_register( 'sendget int scal kind=4 to array kind=4 to image 2 failed.')

      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= co_int_scal_src_k1)) &
        call print_and_register( 'sendget int scal kind=1 to array kind=1 to image 2 failed.')

      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - co_real_scal_src_k8) > tolerance8)) &
        & call print_and_register( 'sendget real kind=8 to array kind=8 to image 2 failed.')

      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - co_real_scal_src_k4) > tolerance4)) &
        & call print_and_register( 'sendget real kind=4 to array kind=4 to image 2 failed.')
    end if

    ! and with kind conversion
    sync all
    if (me == 2) then
      co_int_dst_k4(:)[1] = co_int_scal_src_k1[3]

      co_int_dst_k1(:)[1] = co_int_scal_src_k4[3]
      co_real_dst_k8(:)[1] = co_real_scal_src_k4[3]

      co_real_dst_k4(:)[1] = co_real_scal_src_k8[3]
    end if

    sync all
    if (me == 1) then
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= co_int_scal_src_k4)) &
        call print_and_register( 'sendget int scal kind=1 to array kind=4 to image 2 failed.')

      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= co_int_scal_src_k1)) &
        call print_and_register( 'sendget int scal kind=4 to array kind=1 to image 2 failed.')

      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - co_real_scal_src_k8) > tolerance8)) &
        & call print_and_register( 'sendget real kind=4 to array kind=8 to image 2 failed.')

      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - co_real_scal_src_k4) > tolerance4)) &
        & call print_and_register( 'sendget real kind=8 to array kind=4 to image 2 failed.')
    end if
    ! and with type conversion
    sync all
    if (me == 2) then
      co_int_dst_k4(:)[1] = co_real_scal_src_k4[3]

      co_int_dst_k1(:)[1] = co_real_scal_src_k4[3]

      co_real_dst_k8(:)[1] = co_int_scal_src_k4[3]

      co_real_dst_k4(:)[1] = co_int_scal_src_k4[3]
    end if

    sync all
    if (me == 1) then
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= INT(co_real_scal_src_k4, 4))) &
        & call print_and_register( 'sendget real scal kind=4 to int array kind=4 to image 2 failed.')

      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= INT(co_real_scal_src_k4, 1))) &
        & call print_and_register( 'sendget real scal kind=1 to int array kind=1 to image 2 failed.')

      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - co_int_scal_src_k4) > tolerance4to8)) &
        & call print_and_register( 'sendget int kind=4 to real array kind=8 to image 2 failed.')

      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - co_int_scal_src_k4) > tolerance4)) &
        & call print_and_register( 'sendget int kind=4 to real array kind=4 to image 2 failed.')
    end if

    ! Now with strides
    ! First check send/copy to self
    if (me == 1) then
      co_int_dst_k1 = -1
      co_int_dst_k1(::2)[1] = co_int_src_k1(1:3)[1]
      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= [co_int_src_k1(1), INT(-1, 1), co_int_src_k1(2), INT(-1, 1), co_int_src_k1(3)])) &
        & call print_and_register( 'sendget strided int kind=1 to kind=1 self failed.')

      co_int_dst_k4 = -1
      co_int_dst_k4(::2)[1] = co_int_src_k4(1:3)[1]
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= [co_int_src_k4(1), -1, co_int_src_k4(2), -1, co_int_src_k4(3)])) &
        & call print_and_register( 'sendget strided int kind=4 to kind=4 self failed.')

      co_int_dst_k4 = -2
      co_int_dst_k4(2::2)[1] = co_int_src_k1(4:5)[1]
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= [-2, co_int_src_k4(4), -2, co_int_src_k4(5), -2])) &
        & call print_and_register( 'sendget strided int kind=1 to kind=4 self failed.')

      co_int_dst_k1 = -2
      co_int_dst_k1(2::2)[1] = co_int_src_k4(4:5)[1]
      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= [INT(-2, 1), co_int_src_k1(4), INT(-2, 1), co_int_src_k1(5), INT(-2, 1)])) &
        & call print_and_register( 'sendget strided int kind=4 to kind=1 self failed.')

      co_int_dst_k1(1:5) = co_int_src_k1(5:1:-1)
      co_int_dst_k1(::2)[1] = co_int_dst_k1(3:1:-1)[1]
      print *, co_int_dst_k1
      ! Note, indezes two times reversed!
      if (any(co_int_dst_k1 /= [co_int_src_k1(3), co_int_src_k1(4), co_int_src_k1(4), co_int_src_k1(2), co_int_src_k1(5)])) &
        & call print_and_register( 'sendget strided with temp int kind=1 to kind=1 self failed.')

      co_int_dst_k4(1:5) = co_int_src_k4(5:1:-1)
      co_int_dst_k4(::2)[1] = co_int_dst_k4(3:1:-1)[1]
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= [co_int_src_k4(3), co_int_src_k4(4), co_int_src_k4(4), co_int_src_k4(2), co_int_src_k4(5)])) &
       & call print_and_register( 'sendget strided with temp int kind=4 to kind=4 self failed.')
    else if (me == 2) then ! Do the real copy to self checks on image 2
      co_real_dst_k4 = -1.0
      co_real_dst_k4(::2)[2] = co_real_src_k4(1:3)[2]
      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - [co_real_src_k4(1), -1.0, co_real_src_k4(2), -1.0, co_real_src_k4(3)]) > tolerance4)) &
        & call print_and_register( 'sendget strided real kind=4 to kind=4 self failed.')

      co_real_dst_k8 = -1.0
      co_real_dst_k8(::2)[2] = co_real_src_k8(1:3)[2]
      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - [co_real_src_k8(1), REAL(-1.0, 8), co_real_src_k8(2), REAL(-1.0, 8), co_real_src_k8(3)]) &
        & > tolerance8)) call print_and_register( 'sendget strided real kind=8 to kind=8 self failed.')

      co_real_dst_k8 = -2.0
      co_real_dst_k8(2::2)[2] = co_real_src_k4(4:5)[2]
      print *, co_real_dst_k8(1:5), lbound(co_real_dst_k8, 1)
      if (any(abs(co_real_dst_k8 - [REAL(-2.0, 8), co_real_src_k8(4), REAL(-2.0, 8), co_real_src_k8(5), REAL(-2.0, 8)]) &
        & > tolerance4to8)) call print_and_register( 'sendget strided real kind=4 to kind=8 self failed.')

      co_real_dst_k4 = -2.0
      co_real_dst_k4(2::2)[2] = co_real_src_k8(1:2)[2]
      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - [-2.0, co_real_src_k4(1), -2.0, co_real_src_k4(2), -2.0]) > tolerance4)) &
        & call print_and_register( 'sendget strided real kind=8 to kind=4 self failed.')

      co_real_dst_k4(1:5) = co_real_src_k4(5:1:-1)
      co_real_dst_k4(::2)[2] = co_real_dst_k4(3:1:-1)[2]
      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - [co_real_src_k4(3), co_real_src_k4(4), co_real_src_k4(4), co_real_src_k4(2), &
        & co_real_src_k4(5)]) > tolerance4)) &
        & call print_and_register( 'sendget strided with temp real kind=4 to kind=4 self failed.')

      co_real_dst_k8(1:5) = co_real_src_k8(5:1:-1)
      co_real_dst_k8(::2)[2] = co_real_dst_k8(3:1:-1)[2]
      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - [co_real_src_k8(3), co_real_src_k8(4), co_real_src_k8(4), co_real_src_k8(2), &
        & co_real_src_k8(5)]) > tolerance8)) &
        & call print_and_register( 'sendget strided with temp real kind=8 to kind=8 self failed.')
    end if

    ! Transfer to other image now.
    sync all
    co_int_dst_k4 = -1
    co_int_dst_k1 = INT(-1, 1)
    co_real_dst_k8 = -1.0
    co_real_dst_k4 = REAL(-1.0, 4)
    sync all
    if (me == 2) then
      co_int_dst_k4(::2)[1] = co_int_src_k4(1:3)[3]

      co_int_dst_k1(::2)[1] = co_int_src_k1(1:3)[3]

      co_real_dst_k8(::2)[1] = co_real_src_k8(1:3)[3]

      co_real_dst_k4(::2)[1] = co_real_src_k4(1:3)[3]
    end if

    sync all
    if (me == 1) then
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= [co_int_src_k4(1), -1, co_int_src_k4(2), -1, co_int_src_k4(3)])) &
        & call print_and_register( 'strided sendget int kind=4 to kind=4 from image 3 to image 1 failed.')

      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= [co_int_src_k1(1), INT(-1, 1), co_int_src_k1(2), INT(-1, 1), co_int_src_k1(3)])) &
        & call print_and_register( 'strided sendget int kind=1 to kind=1 from image 3 to image 1 failed.')

      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - [co_real_src_k8(1), REAL(-1.0, 8), co_real_src_k8(2), REAL(-1.0, 8), &
        & co_real_src_k8(3)]) > tolerance8)) &
        & call print_and_register( 'strided sendget real kind=8 to kind=8 from image 3 to image 1 failed.')

      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - [co_real_src_k4(1), REAL(-1.0, 4), co_real_src_k4(2), REAL(-1.0, 4), &
        & co_real_src_k4(3)]) > tolerance4)) &
        & call print_and_register( 'strided sendget real kind=4 to kind=4 from image 3 to image 1 failed.')
    end if

    ! now with strides and kind conversion
    sync all
    co_int_dst_k4 = -1
    co_int_dst_k1 = INT(-1, 1)
    co_real_dst_k8 = -1.0
    co_real_dst_k4 = REAL(-1.0, 4)
    sync all
    if (me == 2) then
      co_int_dst_k4(::2)[1] = co_int_src_k1(3:5)[3]

      co_int_dst_k1(::2)[1] = co_int_src_k4(3:5)[3]

      co_real_dst_k8(::2)[1] = co_real_src_k4(3:5)[3]

      co_real_dst_k4(::2)[1] = co_real_src_k8(3:5)[3]
    end if

    sync all
    if (me == 1) then
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= [co_int_src_k4(3), -1, co_int_src_k4(4), -1, co_int_src_k4(5)])) &
        & call print_and_register( 'strided sendget int kind=1 to kind=4 from image 3 to image 1 failed.')

      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= [co_int_src_k1(3), INT(-1, 1), co_int_src_k1(4), INT(-1, 1), co_int_src_k1(5)])) &
        & call print_and_register( 'strided sendget int kind=4 to kind=1 from image 3 to image 1 failed.')

      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - [co_real_src_k8(3), REAL(-1.0, 8), co_real_src_k8(4), REAL(-1.0, 8), &
        & co_real_src_k8(5)]) > tolerance8)) &
        & call print_and_register( 'strided sendget real kind=4 to kind=8 from image 3 to image 1 failed.')

      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - [co_real_src_k4(3), REAL(-1.0, 4), co_real_src_k4(4), REAL(-1.0, 4), &
        & co_real_src_k4(5)]) > tolerance4)) &
        & call print_and_register( 'strided sendget real kind=8 to kind=4 from image 3 to image 1 failed.')
    end if

    ! now with strides and type conversion
    sync all
    co_int_dst_k4 = -1
    co_int_dst_k1 = INT(-1, 1)
    co_real_dst_k8 = -1.0
    co_real_dst_k4 = REAL(-1.0, 4)
    sync all
    if (me == 2) then
      co_int_dst_k4(::2)[1] = co_real_src_k8(1:3)[3]

      co_int_dst_k1(::2)[1] = co_real_src_k4(1:3)[3]

      co_real_dst_k8(::2)[1] = co_int_src_k4(3:5)[3]

      co_real_dst_k4(::2)[1] = co_int_src_k1(3:5)[3]
    end if

    sync all
    if (me == 1) then
      print *, co_int_dst_k4
      if (any(co_int_dst_k4 /= [5, -1, 4, -1, 3])) &
        & call print_and_register( 'strided sendget real kind=4 to int kind=4 from image 3 to image 1 failed.')

      print *, co_int_dst_k1
      if (any(co_int_dst_k1 /= [INT(5, 1), INT(-1, 1), INT(4, 1), INT(-1, 1), INT(3, 1)])) &
        & call print_and_register( 'strided sendget int real kind=4 to int kind=1 from image 3 to image 1 failed.')

      print *, co_real_dst_k8
      if (any(abs(co_real_dst_k8 - [3.0, -1.0, 2.0, -1.0, 1.0]) > tolerance8)) &
        & call print_and_register( 'strided sendget int kind=4 to real kind=8 from image 3 to image 1 failed.')

      print *, co_real_dst_k4
      if (any(abs(co_real_dst_k4 - [REAL(3, 4), REAL(-1.0, 4), REAL(2, 4), REAL(-1.0, 4), REAL(1, 4)]) > tolerance4)) &
        & call print_and_register( 'strided sendget int kind=4 to real kind=4 from image 3 to image 1 failed.')
    end if

    select case(me)
     case(1)
       if (error_printed) error stop
       sync images([2,3])
       print *, "Test passed."
     case(2)
       if (error_printed) error stop
       sync images(1)
     case(3)
       if (error_printed) error stop
       sync images(1)
    end select

  end associate

contains

  subroutine print_and_register(error_message)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_message
    write(error_unit,*) error_message
    error_printed=.true.
  end subroutine

end program sendget_convert_nums

! vim:ts=2:sts=2:sw=2:
