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
  use iso_fortran_env, only : int8,int32,real32,real64

  implicit none

  real(kind=real32), parameter :: tolerance4 = 1.0e-4_real32
  real(kind=real64), parameter :: tolerance4to8 = 1.0E-4_real64
  real(kind=real64), parameter :: tolerance8 = 1.0E-6_real64

  type t
    integer(kind=int8), allocatable :: int_scal_k1
    integer(kind=int32), allocatable :: int_scal_k4
    real(kind=real32)   , allocatable :: real_scal_k4
    real(kind=real64)   , allocatable :: real_scal_k8
    integer(kind=int8), allocatable, dimension(:) :: int_k1
    integer(kind=int32), allocatable, dimension(:) :: int_k4
    real(kind=real32)   , allocatable, dimension(:) :: real_k4
    real(kind=real64)   , allocatable, dimension(:) :: real_k8
  end type t

  integer(kind=int8)                              :: int_scal_k1
  integer(kind=int32)                              :: int_scal_k4
  real(kind=real32)                                 :: real_scal_k4
  real(kind=real64)                                 :: real_scal_k8

  integer(kind=int8), dimension(1:5)                            :: int_k1
  integer(kind=int32), dimension(1:5)                            :: int_k4
  real(kind=real32)   , dimension(1:5)                            :: real_k4
  real(kind=real64)   , dimension(1:5)                            :: real_k8

  type(t), save, codimension[*] :: obj

  logical :: error_printed=.false.

  associate(me => this_image(), np => num_images())
    if (np < 2) error stop 'Cannot run with less than 2 images.'

    int_scal_k1 = INT(42, kind(int_scal_k1))
    int_scal_k4 = 42
    int_k1 = INT([5, 4, 3, 2, 1], kind(int_scal_k4))
    int_k4 = [5, 4, 3, 2, 1]
    allocate(obj%int_scal_k1, obj%int_scal_k4, obj%int_k1(5), obj%int_k4(5)) ! allocate syncs here

    real_scal_k4 = 37.042
    real_scal_k8 = REAL(37.042, kind(real_scal_k8))
    real_k4 = [ 5.1, 4.2, 3.3, 2.4, 1.5]
    real_k8 = REAL([ 5.1, 4.2, 3.3, 2.4, 1.5], kind(real_k8))
    allocate(obj%real_scal_k4, obj%real_scal_k8, obj%real_k4(1:5), obj%real_k8(1:5)) ! allocate syncs here

    ! Now with strides
    ! First check send/copy to self
    if (me == 1) then

      obj%int_k1(1:5) = int_k1(5:1:-1)
      obj[1]%int_k1(::2) = obj%int_k1(3:1:-1)
      print *, obj%int_k1
      ! Note, indezes two times reversed!
      if (any(obj%int_k1 /= [int_k1(3), int_k1(4), int_k1(4), int_k1(2), int_k1(5)])) &
        & call print_and_register( 'send strided with temp int kind=int8 to kind=1 self failed')

       obj%int_k4(1:5) = int_k4(5:1:-1)
       obj[1]%int_k4(::2) = obj%int_k4(3:1:-1)
       print *, obj%int_k4
       if (any(obj%int_k4 /= [int_k4(3), int_k4(4), int_k4(4), int_k4(2), int_k4(5)])) &
        & call print_and_register( 'send strided with temp int kind=4 to kind=4 self failed')
    else if (me == 2) then ! Do the real copy to self checks on image 2

       obj%real_k4(1:5) = real_k4(5:1:-1)
       obj[2]%real_k4(::2) = obj%real_k4(3:1:-1)
       print *, obj%real_k4
       if (any(abs(obj%real_k4 - [real_k4(3), real_k4(4), real_k4(4), real_k4(2), real_k4(5)]) > tolerance4)) &
         & call print_and_register( 'send strided with temp real kind=4 to kind=4 self failed')

       obj%real_k8(1:5) = real_k8(5:1:-1)
       obj[2]%real_k8(::2) = obj%real_k8(3:1:-1)
       print *, obj%real_k8
       if (any(abs(obj%real_k8 - [real_k8(3), real_k8(4), real_k8(4), real_k8(2), real_k8(5)]) > tolerance8)) &
         & call print_and_register( 'send strided with temp real kind=real64 to kind=8 self failed')
    end if

    select case(me)
      case(1)
         sync images(2) ! wait for image 2 to finish all checks
         if (error_printed) error stop
         sync images(2) ! wait for image 2 to get past its conditional error termination
         print *, "Test passed."
      case(2)
         sync images(1) ! wait for image 1 to finish all checks
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
