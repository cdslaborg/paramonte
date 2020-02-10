!! Thoroughly test get, i.e. foo = bar[N] in all variants
!!
!! Do simple tests for get(). These test comprise
!!
!! FOO = BAR [N]
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

program get_convert_nums

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
    if (np < 2) error stop 'Can not run with less than 2 images.'

    allocate(obj%int_scal_k1, SOURCE=INT(42, 1)) ! allocate syncs here
    allocate(obj%int_scal_k4, SOURCE=42) ! allocate syncs here
    allocate(obj%int_k1(5), SOURCE=INT([ 5, 4, 3, 2, 1], 1)) ! allocate syncs here
    allocate(obj%int_k4(5), SOURCE=[ 5, 4, 3, 2, 1]) ! allocate syncs here

    allocate(obj%real_scal_k4, SOURCE=37.042) ! allocate syncs here
    allocate(obj%real_scal_k8, SOURCE=REAL(37.042, 8)) ! allocate syncs here
    allocate(obj%real_k4(1:5), SOURCE=[ 5.1, 4.2, 3.3, 2.4, 1.5]) ! allocate syncs here
    allocate(obj%real_k8(1:5), SOURCE=REAL([ 5.1, 4.2, 3.3, 2.4, 1.5], 8)) ! allocate syncs here

    ! First check send/copy to self
    if (me == 1) then
      int_scal_k1 = obj[1]%int_scal_k1
      print *, int_scal_k1
      if (obj%int_scal_k1 /= int_scal_k1) call print_and_register('get scalar int kind=1 from kind=1 self failed.')

      int_scal_k4 = obj[1]%int_scal_k4
      print *, int_scal_k4
      if (obj%int_scal_k4 /= int_scal_k4) call print_and_register( 'get scalar int kind=4 to kind=4 self failed.')

      int_scal_k4 = obj[1]%int_scal_k1
      print *, int_scal_k4
      if (obj%int_scal_k4 /= int_scal_k4) call print_and_register( 'get scalar int kind=1 to kind=4 self failed.')

      int_scal_k1 = obj[1]%int_scal_k4
      print *, int_scal_k1
      if (obj%int_scal_k1 /= int_scal_k1) call print_and_register( 'get scalar int kind=4 to kind=1 self failed.')

      int_k1(:) = obj[1]%int_k1(:)
      print *, int_k1
      if (any(obj%int_k1 /= int_k1)) call print_and_register( 'get int kind=1 to kind=1 self failed.')

      int_k4(:) = obj[1]%int_k4(:)
      print *, int_k4
      if (any(obj%int_k4 /= int_k4)) call print_and_register( 'get int kind=4 to kind=4 self failed.')

      int_k4(:) = obj[1]%int_k1(:)
      print *, int_k4
      if (any(obj%int_k4 /= int_k4)) call print_and_register( 'get int kind=1 to kind=4 self failed.')

      int_k1(:) = obj[1]%int_k4(:)
      print *, int_k1
      if (any(obj%int_k1 /= int_k1)) call print_and_register( 'get int kind=4 to kind=1 self failed.')
    else if (me == 2) then ! Do the real copy to self checks on image 2
      real_scal_k4 = obj[2]%real_scal_k4
      print *, real_scal_k4
      if (abs(obj%real_scal_k4 - real_scal_k4) > tolerance4) &
        call print_and_register( 'get scalar real kind=4 to kind=4 self failed.')


      real_scal_k8 = obj[2]%real_scal_k8
      print *, real_scal_k8
      if (abs(obj%real_scal_k8 - real_scal_k8) > tolerance8) &
        call print_and_register( 'get scalar real kind=8 to kind=8 self failed.')


      real_scal_k8 = obj[2]%real_scal_k4
      print *, real_scal_k8
      if (abs(obj%real_scal_k8 - real_scal_k8) > tolerance4to8) &
        call print_and_register( 'get scalar real kind=4 to kind=8 self failed.')


      real_scal_k4 = obj[2]%real_scal_k8
      print *, real_scal_k4
      if (abs(obj%real_scal_k4 - real_scal_k4) > tolerance4) &
        call print_and_register( 'get scalar real kind=8 to kind=4 self failed.')

      real_k4(:) = obj[2]%real_k4(:)
      print *, real_k4
      if (any(abs(obj%real_k4 - real_k4) > tolerance4)) call print_and_register( 'get real kind=4 to kind=4 self failed.')

      real_k8(:) = obj[2]%real_k8(:)
      print *, real_k8
      if (any(abs(obj%real_k8 - real_k8) > tolerance8)) call print_and_register( 'get real kind=8 to kind=8 self failed.')

      real_k8(:) = obj[2]%real_k4(:)
      print *, real_k8
      if (any(abs(obj%real_k8 - real_k8) > tolerance4to8)) call print_and_register( 'get real kind=4 to kind=8 self failed.')

      real_k4(:) = obj[2]%real_k8(:)
      print *, real_k4
      if (any(abs(obj%real_k4 - real_k4) > tolerance4)) call print_and_register( 'get real kind=8 to kind=4 self failed.')
    end if

    sync all
    if (me == 1) then
      int_scal_k1 = obj[2]%int_scal_k1
      print *, int_scal_k1
      if (obj%int_scal_k1 /= int_scal_k1) call print_and_register( 'get scalar int kind=1 to kind=1 to image 2 failed.')

      int_scal_k4 = obj[2]%int_scal_k4
      print *, int_scal_k4
      if (obj%int_scal_k4 /= int_scal_k4) call print_and_register( 'get scalar int kind=4 to kind=4 to image 2 failed.')

      int_k1(:) = obj[2]%int_k1(:)
      print *, int_k1
      if (any(obj%int_k1 /= int_k1)) call print_and_register( 'get int kind=1 to kind=1 to image 2 failed.')

      int_k4(:) = obj[2]%int_k4(:)
      print *, int_k4
      if (any(obj%int_k4 /= int_k4)) call print_and_register( 'get int kind=4 to kind=4 to image 2 failed.')

    else if (me == 2) then
      real_scal_k4 = obj[1]%real_scal_k4
      print *, real_scal_k4
      if (abs(obj%real_scal_k4 - real_scal_k4) > tolerance4) &
        call print_and_register( 'get scalar real kind=4 to kind=4 to image 2 failed.')

      real_scal_k8 = obj[1]%real_scal_k8
      print *, real_scal_k8
      if (abs(obj%real_scal_k8 - real_scal_k8) > tolerance8) &
        call print_and_register( 'get scalar real kind=8 to kind=8 to image 2 failed.')

      real_k4(:) = obj[1]%real_k4(:)
      print *, real_k4
      if (any(abs(obj%real_k4 - real_k4) > tolerance4)) call print_and_register( 'get real kind=4 to kind=4 to image 2 failed.')

      real_k8(:) = obj[1]%real_k8(:)
      print *, real_k8
      if (any(abs(obj%real_k8 - real_k8) > tolerance8)) call print_and_register( 'get real kind=8 to kind=8 to image 2 failed.')
    end if

    sync all
    if (me == 1) then
      int_scal_k4 = obj[2]%int_scal_k1
      print *, int_scal_k4
      if (obj%int_scal_k4 /= int_scal_k4) call print_and_register( 'get scalar int kind=1 to kind=4 to image 2 failed.')

      int_scal_k1 = obj[2]%int_scal_k4
      print *, int_scal_k1
      if (obj%int_scal_k1 /= int_scal_k1) call print_and_register( 'get scalar int kind=4 to kind=1 to image 2 failed.')

      int_k4(:) = obj[2]%int_k1(:)
      print *, int_k4
      if (any(obj%int_k4 /= int_k4)) call print_and_register( 'get int kind=1 to kind=4 to image 2 failed.')

      int_k1(:) = obj[2]%int_k4(:)
      print *, int_k1
      if (any(obj%int_k1 /= int_k1)) call print_and_register( 'get int kind=4 to kind=1 to image 2 failed.')

    elseif (me == 2) then
      real_scal_k8 = obj[1]%real_scal_k4
      print *, real_scal_k8
      if (abs(obj%real_scal_k8 - real_scal_k8) > tolerance4to8) &
        call print_and_register( 'get scalar real kind=4 to kind=8 to image 2 failed.')

      real_scal_k4 = obj[1]%real_scal_k8
      print *, real_scal_k4
      if (abs(obj%real_scal_k4 - real_scal_k4) > tolerance4) &
        call print_and_register( 'get scalar real kind=8 to kind=4 to image 2 failed.')

      real_k8(:) = obj[1]%real_k4(:)
      print *, real_k8
      if (any(abs(obj%real_k8 - real_k8) > tolerance4to8)) &
        call print_and_register( 'get real kind=4 to kind=8 to image 2 failed.')

      real_k4(:) = obj[1]%real_k8(:)
      print *, real_k4
      if (any(abs(obj%real_k4 - real_k4) > tolerance4)) call print_and_register( 'get real kind=8 to kind=4 to image 2 failed.')
    end if

    ! Scalar to array replication
    sync all
    if (me == 1) then
      int_k4(:) = obj[2]%int_scal_k4
      print *, int_k4
      if (any(obj%int_scal_k4 /= int_k4)) call print_and_register( 'get int scal kind=4 to array kind=4 to image 2 failed.')

      int_k1(:) = obj[2]%int_scal_k1
      print *, int_k1
      if (any(obj%int_scal_k1 /= int_k1)) call print_and_register( 'get int scal kind=1 to array kind=1 to image 2 failed.')

    else if (me == 2) then
      real_k8(:) = obj[1]%real_scal_k8
      print *, real_k8
      if (any(abs(obj%real_scal_k8 - real_k8) > tolerance8)) &
        call print_and_register( 'get real kind=8 to array kind=8 to image 2 failed.')

      real_k4(:) = obj[1]%real_scal_k4
      print *, real_k4
      if (any(abs(obj%real_scal_k4 - real_k4) > tolerance4)) &
        call print_and_register( 'get real kind=4 to array kind=4 to image 2 failed.')
    end if

    ! and with kind conversion
    sync all
    if (me == 1) then
      int_k4(:) = obj[2]%int_scal_k1
      print *, int_k4
      if (any(obj%int_scal_k4 /= int_k4)) call print_and_register( 'get int scal kind=1 to array kind=4 to image 2 failed.')

      int_k1(:) = obj[2]%int_scal_k4
      print *, int_k1
      if (any(obj%int_scal_k1 /= int_k1)) call print_and_register( 'get int scal kind=4 to array kind=1 to image 2 failed.')

    else if (me == 2) then
      real_k8(:) = obj[1]%real_scal_k4
      print *, real_k8
      if (any(abs(obj%real_scal_k8 - real_k8) > tolerance8)) &
        call print_and_register( 'get real kind=4 to array kind=8 to image 2 failed.')

      real_k4(:) = obj[1]%real_scal_k8
      print *, real_k4
      if (any(abs(obj%real_scal_k4 - real_k4) > tolerance4)) &
        call print_and_register( 'get real kind=8 to array kind=4 to image 2 failed.')
    end if

    ! and with type conversion
    sync all
    if (me == 1) then
      int_k4(:) = obj[2]%real_scal_k4
      print *, int_k4
      if (any(int_k4 /= INT(obj%real_scal_k4, 4))) &
        call print_and_register( 'get real scal kind=4 to int array kind=4 to image 2 failed.')

      int_k1(:) = obj[2]%real_scal_k4
      print *, int_k1
      if (any(int_k1 /= INT(obj%real_scal_k4, 1))) &
        call print_and_register( 'get real scal kind=1 to int array kind=1 to image 2 failed.')

    else if (me == 2) then
      real_k8(:) = obj[1]%int_scal_k4
      print *, real_k8
      if (any(abs(real_k8 - obj%int_scal_k4) > tolerance4to8)) &
        call print_and_register( 'get int kind=4 to real array kind=8 to image 2 failed.')


      real_k4(:) = obj[1]%int_scal_k4
      print *, real_k4
      if (any(abs(real_k4 - obj%int_scal_k4) > tolerance4)) &
        call print_and_register( 'get int kind=4 to real array kind=4 to image 2 failed.')
    end if

    sync all

    ! Now with strides

    ! Transfer to other image now.
    sync all
    int_k4 = -1
    int_k1 = INT(-1, 1)
    real_k8 = -1.0
    real_k4 = REAL(-1.0, 4)

    sync all
    if (me == 1) then
      int_k4(1:3) = obj[2]%int_k4(::2)
      print *, int_k4
      if (any(int_k4 /= [obj%int_k4(1), obj%int_k4(3), obj%int_k4(5), -1, -1])) &
        & call print_and_register( 'strided get int kind=4 to kind=4 to image 2 failed.')

      int_k1(3:5) = obj[2]%int_k1(::2)
      print *, int_k1
      if (any(int_k1 /= [INT(-1, 1), INT(-1, 1), obj%int_k1(1), obj%int_k1(3), obj%int_k1(5)])) &
        & call print_and_register( 'strided get int kind=1 to kind=1 to image 2 failed.')

      real_k8(1:3) = obj[2]%real_k8(::2)
      print *, real_k8
      if (any(abs(real_k8 - [obj%real_k8(1), obj%real_k8(3), obj%real_k8(5), REAL(-1.0, 8), REAL(-1.0, 8)]) > tolerance8)) &
        & call print_and_register( 'strided get real kind=8 to kind=8 to image 2 failed.')

      real_k4(3:5) = obj[2]%real_k4(::2)
      print *, real_k4
      if (any(abs(real_k4 - [-1.0, -1.0, obj%real_k4(1), obj%real_k4(3), obj%real_k4(5)]) > tolerance4)) &
        & call print_and_register( 'strided get real kind=4 to kind=4 to image 2 failed.')
    end if

    ! now with strides and kind conversion
    sync all
    int_k4 = -1
    int_k1 = -1
    real_k4 = -1.0
    real_k8 = -1.0
    obj%int_k4 = [105, 104, 103, 102, 101]
    obj%int_k1 = INT([15, 14, 13, 12, 11], 1)
    obj%real_k8 = [5.1, 4.2, 3.3, 2.4, 1.5]
    obj%real_k4 = REAL([-5.1, -4.2, -3.3, -2.4, -1.5], 4)
    sync all
    if (me == 1) then
      int_k4(1:3) = obj[2]%int_k1(::2)
      print *, int_k4
      if (any(int_k4 /= [15, 13, 11, -1, -1])) call print_and_register( 'strided get int kind=1 to kind=4 to image 2 failed.')

      int_k1(1:3) = obj[2]%int_k4(::2)
      print *, int_k1
      if (any(int_k1 /= INT([105, 103, 101, -1, -1], 1))) &
        & call print_and_register( 'strided get int kind=4 to kind=1 to image 2 failed.')

      real_k8(1:3) = obj[2]%real_k4(::2)
      print *, real_k8
      if (any(abs(real_k8 - [-5.1, -3.3, -1.5, -1.0, -1.0]) > tolerance8)) &
        & call print_and_register( 'strided get real kind=4 to kind=8 to image 2 failed.')

      real_k4(1:3) = obj[2]%real_k8(::2)
      print *, real_k4
      if (any(abs(real_k4 - REAL([5.1, 3.3, 1.5, -1.0, -1.0], 4)) > tolerance4)) &
        & call print_and_register( 'strided get real kind=8 to kind=4 to image 2 failed.')

    else if (me == 2) then
      ! now with strides and type conversion
      int_k4(1:3) = obj[1]%real_k8(::2)
      print *, int_k4
      if (any(int_k4 /= [5, 3, 1, -1, -1])) call print_and_register( 'strided get real kind=4 to int kind=4 to image 2 failed.')

      int_k1(1:3) = obj[1]%real_k4(::2)
      print *, int_k1
      if (any(int_k1 /= INT([-5, -3, -1, -1, -1], 1))) &
        & call print_and_register( 'strided get real kind=4 to int kind=1 to image 2 failed.')

      real_k8(1:3) = obj[1]%int_k4(::2)
      print *, real_k8
      if (any(abs(real_k8 - [105.0, 103.0, 101.0, -1.0, -1.0]) > tolerance8)) &
        & call print_and_register( 'strided get int kind=4 to real kind=8 to image 2 failed.')

      real_k4(1:3) = obj[1]%int_k1(::2)
      print *, real_k4
      if (any(abs(real_k4 - [15.0, 13.0, 11.0, -1.0, -1.0]) > tolerance4)) &
        & call print_and_register( 'strided get int kind=1 to real kind=4 to image 2 failed.')
    end if

    if (error_printed) error stop
    if (me==2) sync images(1)
    if (me==1) then
      sync images(2)
      print *, "Test passed."
    end if
  end associate


contains

  subroutine print_and_register(error_message)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: error_message
    write(error_unit,*) error_message
    error_printed=.true.
  end subroutine

end program get_convert_nums

! vim:ts=2:sts=2:sw=2:
