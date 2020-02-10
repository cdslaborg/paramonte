!! Thoroughly test send, i.e. foo [N] = bar in all variants
!!
!! Do simple tests for send(). These test comprise
!!
!! FOO [N] = BAR
!!
!! where
!!
!!  FOO                BAR                 images
!! character(len=20) character(len=10)   N == me
!!  kind == 1         kind == 1
!!  kind == 1         kind == 4
!!  kind == 4         kind == 1
!!  kind == 4         kind == 4
!!
!! character         character
!!   (1:4, len == 5)  (1:4, len == 5)
!!  kind == 1         kind == 1
!!  kind == 1         kind == 4
!!  kind == 4         kind == 1
!!  kind == 4         kind == 4
!!
!! all of the above but for              N != me
!!
!! And may be some other, I've forgotten.
!!
!! Author: Andre Vehreschild, 2017



program send_convert_char_array

  implicit none

  character(kind=1, len=:), allocatable, codimension[:] :: co_str_k1_scal
  character(kind=1, len=:), allocatable :: str_k1_scal
  character(kind=4, len=:), allocatable, codimension[:] :: co_str_k4_scal
  character(kind=4, len=:), allocatable :: str_k4_scal

  character(kind=1, len=:), allocatable, codimension[:] :: co_str_k1_arr(:)
  character(kind=1, len=:), allocatable :: str_k1_arr(:)
  character(kind=4, len=:), allocatable, codimension[:] :: co_str_k4_arr(:)
  character(kind=4, len=:), allocatable :: str_k4_arr(:)

  logical ::  error_printed=.false.

  associate(me => this_image(), np => num_images())
    if (np < 2) error stop 'Can not run with less than 2 images.'

    allocate(str_k1_scal, SOURCE='abcdefghij')
    allocate(str_k4_scal, SOURCE=4_'abcdefghij')
    allocate(character(len=20)::co_str_k1_scal[*]) ! allocate syncs here
    allocate(character(kind=4, len=20)::co_str_k4_scal[*]) ! allocate syncs here

    allocate(str_k1_arr(1:4), SOURCE=['abc', 'EFG', 'klm', 'NOP'])
    allocate(str_k4_arr(1:4), SOURCE=[4_'abc', 4_'EFG', 4_'klm', 4_'NOP'])
    allocate(character(len=5)::co_str_k1_arr(4)[*])
    allocate(character(kind=4, len=5)::co_str_k4_arr(4)[*])

    ! First check send/copy to self
    if (me == 1) then
      co_str_k1_scal[1] = str_k1_scal
      print *, '#' // co_str_k1_scal // '#, len:', len(co_str_k1_scal)
      if (co_str_k1_scal /= str_k1_scal // '          ') call print_and_register( 'send scalar kind=1 to kind=1 self failed.')

      co_str_k4_scal[1] = str_k4_scal
      print *, 4_'#' // co_str_k4_scal // 4_'#, len:', len(co_str_k4_scal)
      if (co_str_k4_scal /= str_k4_scal // 4_'          ') call print_and_register( 'send scalar kind=4 to kind=4 self failed.')

      co_str_k4_scal[1] = str_k1_scal
      print *, 4_'#' // co_str_k4_scal // 4_'#, len:', len(co_str_k4_scal)
      if (co_str_k4_scal /= str_k4_scal // 4_'          ') call print_and_register( 'send scalar kind=1 to kind=4 self failed.')

      co_str_k1_scal[1] = str_k4_scal
      print *, '#' // co_str_k1_scal // '#, len:', len(co_str_k1_scal)
      if (co_str_k1_scal /= str_k1_scal // '          ') call print_and_register( 'send scalar kind=4 to kind=1 self failed.')
    end if

    ! Do the same for arrays but on image 2
    if (me == 2) then
      co_str_k1_arr(:)[2] = str_k1_arr
      print *, '#' // co_str_k1_arr(:) // '#, len:', len(co_str_k1_arr(1))
      if (any(co_str_k1_arr /= ['abc  ', 'EFG  ', 'klm  ', 'NOP  '])) &
        call print_and_register( 'send array kind=1 to kind=1 self failed.')

      print *, str_k4_arr
      co_str_k4_arr(:)[2] = [4_'abc', 4_'EFG', 4_'klm', 4_'NOP']! str_k4_arr
      print *, 4_'#' // co_str_k4_arr(:) // 4_'#, len:', len(co_str_k4_arr(1))
      if (any(co_str_k4_arr /= [4_'abc  ', 4_'EFG  ', 4_'klm  ', 4_'NOP  '])) &
        call print_and_register( 'send array kind=4 to kind=4 self failed.')

      co_str_k4_arr(:)[2] = str_k1_arr
      print *, 4_'#' // co_str_k4_arr(:) // 4_'#, len:', len(co_str_k4_arr(1))
      if (any(co_str_k4_arr /= [ 4_'abc  ', 4_'EFG  ', 4_'klm  ', 4_'NOP  '])) &
        call print_and_register( 'send array kind=1 to kind=4 self failed.')

      co_str_k1_arr(:)[2] = str_k4_arr
      print *, '#' // co_str_k1_arr(:) // '#, len:', len(co_str_k1_arr(1))
      if (any(co_str_k1_arr /= ['abc  ', 'EFG  ', 'klm  ', 'NOP  '])) &
        call print_and_register( 'send array kind=4 to kind=1 self failed.')
    end if

    sync all
    if (me == 1) then
      co_str_k1_scal[2] = str_k1_scal

      co_str_k4_scal[2] = str_k4_scal
    end if

    sync all
    if (me == 2) then
      print *, '#' // co_str_k1_scal // '#, len:', len(co_str_k1_scal)
      if (co_str_k1_scal /= str_k1_scal // '          ') call print_and_register( 'send kind=1 to kind=1 image 2 failed.')

      print *, 4_'#' // co_str_k4_scal // 4_'#, len:', len(co_str_k4_scal)
      if (co_str_k4_scal /= str_k4_scal // 4_'          ') call print_and_register( 'send kind=4 to kind=4 image 2 failed.')
    end if

    sync all
    if (me == 1) then
      co_str_k4_scal[2] = str_k1_scal

      co_str_k1_scal[2] = str_k4_scal
    end if

    sync all
    if (me == 2) then
      print *, 4_'#' // co_str_k4_scal // 4_'#, len:', len(co_str_k4_scal)
      if (co_str_k4_scal /= str_k4_scal // 4_'          ') call print_and_register( 'send kind=1 to kind=4 to image 2 failed.')

      print *, '#' // co_str_k1_scal // '#, len:', len(co_str_k1_scal)
      if (co_str_k1_scal /= str_k1_scal // '          ') call print_and_register( 'send kind=4 to kind=1 to image 2 failed.')
    end if

    co_str_k1_arr(:) = '#####'
    co_str_k4_arr(:) = 4_'#####'

    sync all

    if (me == 1) then
      co_str_k1_arr(::2)[2] = 'foo'
      co_str_k4_arr(::2)[2] = ['bar', 'baz']
    end if

    sync all
    if (me == 2) then
      print *, co_str_k1_arr
      if (any(co_str_k1_arr /= ['foo  ', '#####', 'foo  ', '#####'])) &
        & call print_and_register( "strided send char arr kind 1 to kind 1 failed.")
      print *, co_str_k4_arr
      if (any(co_str_k4_arr /= [4_'bar  ', 4_'#####', 4_'baz  ', 4_'#####'] )) &
        & call print_and_register( "strided send char arr kind 1 to kind 4 failed.")
    end if

    select case(me)
      case(1)
        if (error_printed) error stop
        sync images(2)
        print *, 'Test passed.'
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

end program send_convert_char_array

! vim:ts=2:sts=2:sw=2:
