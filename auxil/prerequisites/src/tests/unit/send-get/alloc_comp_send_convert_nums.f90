!! Thoroughly test send, i.e. foo[N].comp = bar in all variants
!!
!! Do simple tests for send(). These test comprise
!!
!! FOO[N].COMP = BAR
!!
!! where
!!
!!  FOO                BAR              images
!! scalar            scalar            N == me
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
!! all of the above but for            N != me
!!
!! And may be some other, I've forgotten.
!!
!! Author: Andre Vehreschild, 2017

program alloc_comp_send_convert_nums

  implicit none

  real(kind=4), parameter :: tolerance4 = 1.0e-4
  real(kind=8), parameter :: tolerance4to8 = 1.0E-4
  real(kind=8), parameter :: tolerance8 = 1.0E-6

  type t
    integer(kind=1), allocatable :: int_scal_k1
    integer(kind=4), allocatable :: int_scal_k4
    real(kind=4)   , allocatable :: real_scal_k4
    real(kind=8)   , allocatable :: real_scal_k8
    integer(kind=1), allocatable, dimension(:) :: int_k1
    integer(kind=4), allocatable, dimension(:) :: int_k4
    real(kind=4)   , allocatable, dimension(:) :: real_k4
    real(kind=8)   , allocatable, dimension(:) :: real_k8
  end type t

  integer(kind=1)                              :: int_scal_k1
  integer(kind=4)                              :: int_scal_k4
  real(kind=4)                                 :: real_scal_k4
  real(kind=8)                                 :: real_scal_k8

  integer(kind=1), dimension(1:5)                            :: int_k1
  integer(kind=4), dimension(1:5)                            :: int_k4
  real(kind=4)   , dimension(1:5)                            :: real_k4
  real(kind=8)   , dimension(1:5)                            :: real_k8

  type(t), save, codimension[*] :: obj

  logical :: error_printed=.false.

  associate(me => this_image(), np => num_images())
    if (np < 2) error stop 'Cannot run with less than 2 images.'

    int_scal_k1 = INT(42, 1)
    int_scal_k4 = 42
    int_k1 = INT([5, 4, 3, 2, 1], 1)
    int_k4 = [5, 4, 3, 2, 1]
    allocate(obj%int_scal_k1, obj%int_scal_k4, obj%int_k1(5), obj%int_k4(5)) ! allocate syncs here

    real_scal_k4 = 37.042
    real_scal_k8 = REAL(37.042, 8)
    real_k4 = [ 5.1, 4.2, 3.3, 2.4, 1.5]
    real_k8 = REAL([ 5.1, 4.2, 3.3, 2.4, 1.5], 8)
    allocate(obj%real_scal_k4, obj%real_scal_k8, obj%real_k4(1:5), obj%real_k8(1:5)) ! allocate syncs here
    ! First check send/copy to self
    if (me == 1) then
      obj[1]%int_scal_k1 = int_scal_k1
      print *, obj%int_scal_k1
      if (obj%int_scal_k1 /= int_scal_k1) error stop 'send scalar int kind=1 to kind=1 self failed.'

      obj[1]%int_scal_k4 = int_scal_k4
      print *, obj%int_scal_k4
      if (obj%int_scal_k4 /= int_scal_k4) call print_and_register( 'send scalar int kind=4 to kind=4 self failed.')

      obj[1]%int_scal_k4 = int_scal_k1
      print *, obj%int_scal_k4
      if (obj%int_scal_k4 /= int_scal_k4) call print_and_register( 'send scalar int kind=1 to kind=4 self failed.')

      obj[1]%int_scal_k1 = int_scal_k4
      print *, obj%int_scal_k1
      if (obj%int_scal_k1 /= int_scal_k1) call print_and_register( 'send scalar int kind=4 to kind=1 self failed.')

      obj[1]%int_k1(:) = int_k1
      print *, obj%int_k1
      if (any(obj%int_k1 /= int_k1)) call print_and_register( 'send int kind=1 to kind=1 self failed.')

      obj[1]%int_k4(:) = int_k4
      print *, obj%int_k4
      if (any(obj%int_k4 /= int_k4)) call print_and_register( 'send int kind=4 to kind=4 self failed.')

      obj[1]%int_k4(:) = int_k1
      print *, obj%int_k4
      if (any(obj%int_k4 /= int_k4)) call print_and_register( 'send int kind=1 to kind=4 self failed.')

      obj[1]%int_k1(:) = int_k4
      print *, obj%int_k1
      if (any(obj%int_k1 /= int_k1)) call print_and_register( 'send int kind=4 to kind=1 self failed.')
    else if (me == 2) then ! Do the real copy to self checks on image 2
      obj[2]%real_scal_k4 = real_scal_k4
      print *, obj%real_scal_k4
      if (abs(obj%real_scal_k4 - real_scal_k4) > tolerance4) &
        call print_and_register( 'send scalar real kind=4 to kind=4 self failed.')

      obj[2]%real_scal_k8 = real_scal_k8
      print *, obj%real_scal_k8
      if (abs(obj%real_scal_k8 - real_scal_k8) > tolerance8) &
        call print_and_register( 'send scalar real kind=8 to kind=8 self failed.')

      obj[2]%real_scal_k8 = real_scal_k4
      print *, obj%real_scal_k8
      if (abs(obj%real_scal_k8 - real_scal_k8) > tolerance4to8) &
        call print_and_register( 'send scalar real kind=4 to kind=8 self failed.')

      obj[2]%real_scal_k4 = real_scal_k8
      print *, obj%real_scal_k4
      if (abs(obj%real_scal_k4 - real_scal_k4) > tolerance4) &
        call print_and_register( 'send scalar real kind=8 to kind=4 self failed.')

      obj[2]%real_k4(:) = real_k4
      print *, obj%real_k4
      if (any(abs(obj%real_k4 - real_k4) > tolerance4)) call print_and_register( 'send real kind=4 to kind=4 self failed.')

      obj[2]%real_k8(:) = real_k8
      print *, obj%real_k8
      if (any(abs(obj%real_k8 - real_k8) > tolerance8)) call print_and_register( 'send real kind=8 to kind=8 self failed.')

      obj[2]%real_k8(:) = real_k4
      print *, obj%real_k8
      if (any(abs(obj%real_k8 - real_k8) > tolerance4to8)) call print_and_register( 'send real kind=4 to kind=8 self failed.')

      obj[2]%real_k4(:) = real_k8
      print *, obj%real_k4
      if (any(abs(obj%real_k4 - real_k4) > tolerance4)) call print_and_register( 'send real kind=8 to kind=4 self failed.')
    end if

    sync all
    if (me == 1) then
      obj[2]%int_scal_k1 = int_scal_k1

      obj[2]%int_scal_k4 = int_scal_k4

      obj[2]%int_k1(:) = int_k1

      obj[2]%int_k4(:) = int_k4

      obj[2]%real_scal_k4 = real_scal_k4

      obj[2]%real_scal_k8 = real_scal_k8

      obj[2]%real_k4(:) = real_k4

      obj[2]%real_k8(:) = real_k8
    end if

    sync all
    if (me == 2) then
      print *, obj%int_scal_k1
      if (obj%int_scal_k1 /= int_scal_k1) call print_and_register( 'send scalar int kind=1 to kind=1 to image 2 failed.')

      print *, obj%int_scal_k4
      if (obj%int_scal_k4 /= int_scal_k4) call print_and_register( 'send scalar int kind=4 to kind=4 to image 2 failed.')

      print *, obj%int_k1
      if (any(obj%int_k1 /= int_k1)) call print_and_register( 'send int kind=1 to kind=1 to image 2 failed.')

      print *, obj%int_k4
      if (any(obj%int_k4 /= int_k4)) call print_and_register( 'send int kind=4 to kind=4 to image 2 failed.')

      print *, obj%real_scal_k4
      if (abs(obj%real_scal_k4 - real_scal_k4) > tolerance4) &
        call print_and_register( 'send scalar real kind=4 to kind=4 to image 2 failed.')

      print *, obj%real_scal_k8
      if (abs(obj%real_scal_k8 - real_scal_k8) > tolerance8) &
        call print_and_register( 'send scalar real kind=8 to kind=8 to image 2 failed.')

      print *, obj%real_k4
      if (any(abs(obj%real_k4 - real_k4) > tolerance4)) &
        call print_and_register( 'send real kind=4 to kind=4 to image 2 failed.')

      print *, obj%real_k8
      if (any(abs(obj%real_k8 - real_k8) > tolerance8)) &
        call print_and_register( 'send real kind=8 to kind=8 to image 2 failed.')
    end if

    sync all
    if (me == 1) then
      obj[2]%int_scal_k4 = int_scal_k1

      obj[2]%int_scal_k1 = int_scal_k4

      obj[2]%int_k4(:) = int_k1

      obj[2]%int_k1(:) = int_k4

      obj[2]%real_scal_k8 = real_scal_k4

      obj[2]%real_scal_k4 = real_scal_k8

      obj[2]%real_k8(:) = real_k4

      obj[2]%real_k4(:) = real_k8
    end if

    sync all
    if (me == 2) then
      print *, obj%int_scal_k4
      if (obj%int_scal_k4 /= int_scal_k4) call print_and_register( 'send scalar int kind=1 to kind=4 to image 2 failed.')

      print *, obj%int_scal_k1
      if (obj%int_scal_k1 /= int_scal_k1) call print_and_register( 'send scalar int kind=4 to kind=1 to image 2 failed.')

      print *, obj%int_k4
      if (any(obj%int_k4 /= int_k4)) call print_and_register( 'send int kind=1 to kind=4 to image 2 failed.')

      print *, obj%int_k1
      if (any(obj%int_k1 /= int_k1)) call print_and_register( 'send int kind=4 to kind=1 to image 2 failed.')

      print *, obj%real_scal_k8
      if (abs(obj%real_scal_k8 - real_scal_k8) > tolerance4to8) &
        call print_and_register( 'send scalar real kind=4 to kind=8 to image 2 failed.')

      print *, obj%real_scal_k4
      if (abs(obj%real_scal_k4 - real_scal_k4) > tolerance4) &
        call print_and_register( 'send scalar real kind=8 to kind=4 to image 2 failed')

      print *, obj%real_k8
      if (any(abs(obj%real_k8 - real_k8) > tolerance4to8)) &
        call print_and_register( 'send real kind=4 to kind=8 to image 2 failed')

      print *, obj%real_k4
      if (any(abs(obj%real_k4 - real_k4) > tolerance4)) call print_and_register( 'send real kind=8 to kind=4 to image 2 failed')
    end if

    ! Scalar to array replication
    sync all
    if (me == 1) then
      obj[2]%int_k4(:) = int_scal_k4

      obj[2]%int_k1(:) = int_scal_k1

      obj[2]%real_k8(:) = real_scal_k8

      obj[2]%real_k4(:) = real_scal_k4
    end if

    sync all
    if (me == 2) then
      print *, obj%int_k4
      if (any(obj%int_k4 /= int_scal_k4)) call print_and_register( 'send int scal kind=4 to array kind=4 to image 2 failed')

      print *, obj%int_k1
      if (any(obj%int_k1 /= int_scal_k1)) call print_and_register( 'send int scal kind=1 to array kind=1 to image 2 failed')

      print *, obj%real_k8
      if (any(abs(obj%real_k8 - real_scal_k8) > tolerance8)) &
        call print_and_register( 'send real kind=8 to array kind=8 to image 2 failed')

      print *, obj%real_k4
      if (any(abs(obj%real_k4 - real_scal_k4) > tolerance4)) &
        call print_and_register( 'send real kind=4 to array kind=4 to image 2 failed')
    end if

    ! and with kind conversion
    sync all
    if (me == 1) then
      obj[2]%int_k4(:) = int_scal_k1

      obj[2]%int_k1(:) = int_scal_k4

      obj[2]%real_k8(:) = real_scal_k4

      obj[2]%real_k4(:) = real_scal_k8
    end if

    sync all
    if (me == 2) then
      print *, obj%int_k4
      if (any(obj%int_k4 /= int_scal_k4)) call print_and_register( 'send int scal kind=1 to array kind=4 to image 2 failed')

      print *, obj%int_k1
      if (any(obj%int_k1 /= int_scal_k1)) call print_and_register( 'send int scal kind=4 to array kind=1 to image 2 failed')

      print *, obj%real_k8
      if (any(abs(obj%real_k8 - real_scal_k8) > tolerance8)) &
        call print_and_register( 'send real kind=4 to array kind=8 to image 2 failed')

      print *, obj%real_k4
      if (any(abs(obj%real_k4 - real_scal_k4) > tolerance4)) &
        call print_and_register( 'send real kind=8 to array kind=4 to image 2 failed')
    end if
    ! and with type conversion
    sync all
    if (me == 1) then
      obj[2]%int_k4(:) = real_scal_k4

      obj[2]%int_k1(:) = real_scal_k4

      obj[2]%real_k8(:) = int_scal_k4

      obj[2]%real_k4(:) = int_scal_k4
    end if

    sync all
    if (me == 2) then
      print *, obj%int_k4
      if (any(obj%int_k4 /= INT(real_scal_k4, 4))) &
        call print_and_register( 'send real scal kind=4 to int array kind=4 to image 2 failed')

      print *, obj%int_k1
      if (any(obj%int_k1 /= INT(real_scal_k4, 1))) &
        call print_and_register( 'send real scal kind=1 to int array kind=1 to image 2 failed')

      print *, obj%real_k8
      if (any(abs(obj%real_k8 - int_scal_k4) > tolerance4to8)) &
        call print_and_register( 'send int kind=4 to real array kind=8 to image 2 failed')

      print *, obj%real_k4
      if (any(abs(obj%real_k4 - int_scal_k4) > tolerance4)) &
        call print_and_register( 'send int kind=4 to real array kind=4 to image 2 failed')
    end if

    ! Now with strides
    ! First check send/copy to self
    if (me == 1) then
      obj%int_k1 = -1
      obj[1]%int_k1(::2) = int_k1(1:3)
      print *, obj%int_k1
      if (any(obj%int_k1 /= [int_k1(1), INT(-1, 1), int_k1(2), INT(-1, 1), int_k1(3)])) &
        & call print_and_register( 'send strided int kind=1 to kind=1 self failed')

      obj%int_k4 = -1
      obj[1]%int_k4(::2) = int_k4
      print *, obj%int_k4
      if (any(obj%int_k4 /= [int_k4(1), -1, int_k4(2), -1, int_k4(3)])) &
        call print_and_register( 'send strided int kind=4 to kind=4 self failed')

      obj%int_k4 = -2
      obj[1]%int_k4(2::2) = int_k1(4:5)
      print *, obj%int_k4
      if (any(obj%int_k4 /= [-2, int_k4(4), -2, int_k4(5), -2])) &
        call print_and_register( 'send strided int kind=1 to kind=4 self failed')

      obj%int_k1 = -2
      obj[1]%int_k1(2::2) = int_k4(4:5)
      print *, obj%int_k1
      if (any(obj%int_k1 /= [INT(-2, 1), int_k1(4), INT(-2, 1), int_k1(5), INT(-2, 1)])) &
        & call print_and_register( 'send strided int kind=4 to kind=1 self failed')

    else if (me == 2) then ! Do the real copy to self checks on image 2
      obj%real_k4 = -1.0
      obj[2]%real_k4(::2) = real_k4(1:3)
      print *, obj%real_k4
      if (any(abs(obj%real_k4 - [real_k4(1), -1.0, real_k4(2), -1.0, real_k4(3)]) > tolerance4)) &
        & call print_and_register( 'send strided real kind=4 to kind=4 self failed')

      obj%real_k8 = -1.0
      obj[2]%real_k8(::2) = real_k8(1:3)
      print *, obj%real_k8
      if (any(abs(obj%real_k8 - [real_k8(1), REAL(-1.0, 8), real_k8(2), REAL(-1.0, 8), real_k8(3)]) > tolerance8)) &
        & call print_and_register( 'send strided real kind=8 to kind=8 self failed')

      obj%real_k8 = -2.0
      obj[2]%real_k8(2::2) = real_k4(4:5)
      print *, obj%real_k8(1:5), lbound(real_k8, 1)
      if (any(abs(obj%real_k8 - [REAL(-2.0, 8), real_k8(4), REAL(-2.0, 8), real_k8(5), REAL(-2.0, 8)]) > tolerance4to8)) &
        & call print_and_register( 'send strided real kind=4 to kind=8 self failed')

      obj%real_k4 = -2.0
      obj[2]%real_k4(2::2) = real_k8(1:2)
      print *, obj%real_k4
      if (any(abs(obj%real_k4 - [-2.0, real_k4(1), -2.0, real_k4(2), -2.0]) > tolerance4)) &
        & call print_and_register( 'send strided real kind=8 to kind=4 self failed')

    end if

    ! Transfer to other image now.
    sync all
    obj%int_k4 = -1
    obj%int_k1 = INT(-1, 1)
    obj%real_k8 = -1.0
    obj%real_k4 = REAL(-1.0, 4)
    sync all
    if (me == 1) then
      obj[2]%int_k4(::2) = [ 15, 13, 11]

      obj[2]%int_k1(::2) = [INT(-15, 1), INT(-13, 1), INT(-11, 1)]

      obj[2]%real_k8(::2) = [REAL(1.3, 8), REAL(1.5, 8), REAL(1.7, 8)]

      obj[2]%real_k4(::2) = [REAL(1.3, 4), REAL(1.5, 4), REAL(1.7, 4)]
    end if

    sync all
    if (me == 2) then
      print *, obj%int_k4
      if (any(obj%int_k4 /= [15, -1, 13, -1, 11])) &
        call print_and_register( 'strided send int kind=4 to kind=4 to image 2 failed')

      print *, obj%int_k1
      if (any(obj%int_k1 /= [INT(-15, 1), INT(-1, 1), INT(-13, 1), INT(-1, 1), INT(-11, 1)])) &
        call print_and_register( 'strided send int kind=1 to kind=1 to image 2 failed')

      print *, obj%real_k8
      if (any(abs(obj%real_k8 - [1.3, -1.0, 1.5, -1.0, 1.7]) > tolerance8)) &
        call print_and_register( 'strided send real kind=8 to kind=8 to image 2 failed')

      print *, obj%real_k4
      if (any(abs(obj%real_k4 - [REAL(1.3, 4), REAL(-1.0, 4), REAL(1.5, 4), REAL(-1.0, 4), REAL(1.7, 4)]) > tolerance4)) &
        call print_and_register( 'strided send real kind=4 to kind=4 to image 2 failed')
    end if

    ! now with strides and kind conversion
    sync all
    obj%int_k4 = -1
    obj%int_k1 = INT(-1, 1)
    obj%real_k8 = -1.0
    obj%real_k4 = REAL(-1.0, 4)
    sync all
    if (me == 1) then
      obj[2]%int_k4(::2) = [INT(15, 1), INT(13, 1), INT(11, 1)]

      obj[2]%int_k1(::2) = [-15, -13, -11]

      obj[2]%real_k8(::2) = [REAL(1.3, 4), REAL(1.5, 4), REAL(1.7, 4)]

      obj[2]%real_k4(::2) = [REAL(1.3, 8), REAL(1.5, 8), REAL(1.7, 8)]
    end if

    sync all
    if (me == 2) then
      print *, obj%int_k4
      if (any(obj%int_k4 /= [15, -1, 13, -1, 11])) call print_and_register( 'strided send int kind=1 to kind=4 to image 2 failed')

      print *, obj%int_k1
      if (any(obj%int_k1 /= [INT(-15, 1), INT(-1, 1), INT(-13, 1), INT(-1, 1), INT(-11, 1)])) &
        & call print_and_register( 'strided send int kind=4 to kind=1 to image 2 failed')

      print *, obj%real_k8
      if (any(abs(obj%real_k8 - [1.3, -1.0, 1.5, -1.0, 1.7]) > tolerance8)) &
        & call print_and_register( 'strided send real kind=4 to kind=8 to image 2 failed')

      print *, obj%real_k4
      if (any(abs(obj%real_k4 - [REAL(1.3, 4), REAL(-1.0, 4), REAL(1.5, 4), REAL(-1.0, 4), REAL(1.7, 4)]) > tolerance4)) &
        & call print_and_register( 'strided send real kind=8 to kind=4 to image 2 failed')
    end if

    ! now with strides and type conversion
    sync all
    obj%int_k4 = -1
    obj%int_k1 = INT(-1, 1)
    obj%real_k8 = -1.0
    obj%real_k4 = REAL(-1.0, 4)
    sync all
    if (me == 1) then
      obj[2]%int_k4(::2) = [15.0, 13.0, 11.0]

      obj[2]%int_k1(::2) = [-15.0, -13.0, -11.0]

      obj[2]%real_k8(::2) = [13, 15, 17]

      obj[2]%real_k4(::2) = [23, 25, 27]
    end if

    sync all
    if (me == 2) then
      print *, obj%int_k4
      if (any(obj%int_k4 /= [15, -1, 13, -1, 11])) &
        call print_and_register( 'strided send real kind=4 to int kind=4 to image 2 failed')

      print *, obj%int_k1
      if (any(obj%int_k1 /= [INT(-15, 1), INT(-1, 1), INT(-13, 1), INT(-1, 1), INT(-11, 1)])) &
          call print_and_register( 'strided send real kind=4 to int kind=1 to image 2 failed')

      print *, obj%real_k8
      if (any(abs(obj%real_k8 - [13.0, -1.0, 15.0, -1.0, 17.0]) > tolerance8)) &
          call print_and_register( 'strided send int kind=4 to real kind=8 to image 2 failed')

      print *, obj%real_k4
      if (any(abs(obj%real_k4 - [REAL(23, 4), REAL(-1.0, 4), REAL(25, 4), REAL(-1.0, 4), REAL(27, 4)]) > tolerance4)) &
          call print_and_register( 'strided send int kind=4 to real kind=4 to image 2 failed')
    end if

    select case(me)
      case(1)
         if (error_printed) error stop
         sync images(2)
         print *, "Test passed."
      case(2)
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

end program alloc_comp_send_convert_nums

! vim:ts=2:sts=2:sw=2:
