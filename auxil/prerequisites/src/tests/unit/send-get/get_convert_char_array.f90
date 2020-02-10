!! Thoroughly test get, i.e. foo = bar [N] in all variants
!!
!! Do simple tests for get(). These test comprise
!!
!! FOO = BAR [N]
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



program get_convert_char_array

  implicit none

  character(kind=1, len=10), codimension[*] :: co_str_k1_scal
  character(kind=1, len=20)                 :: str_k1_scal
  character(kind=4, len=10), codimension[*] :: co_str_k4_scal
  character(kind=4, len=20)                 :: str_k4_scal

  character(kind=1, len=5), codimension[*] :: co_str_k1_arr(1:4)
  character(kind=1, len=5)                 :: str_k1_arr(1:4)
  character(kind=4, len=5), codimension[*] :: co_str_k4_arr(1:4)
  character(kind=4, len=5)                 :: str_k4_arr(1:4)
  logical :: error_printed=.false.

  associate(me => this_image(), np => num_images())
    if (np < 2) error stop 'Can not run with less than 2 images.'

    co_str_k1_scal = 'abcdefghij'
    co_str_k4_scal = 4_'abcdefghij'

    co_str_k1_arr(:) = ['abc', 'EFG', 'klm', 'NOP']
    co_str_k4_arr(:) = [4_'abc', 4_'EFG', 4_'klm', 4_'NOP']

    ! First check get/copy to self
    if (me == 1) then
      str_k1_scal = co_str_k1_scal[1]
      print *, '#' // str_k1_scal // '#, len:', len(str_k1_scal)
      if (co_str_k1_scal /= str_k1_scal // '          ') call print_and_register( 'get scalar kind=1 to kind=1 self failed.')

      str_k4_scal = co_str_k4_scal[1]
      print *, 4_'#' // str_k4_scal // 4_'#, len:', len(str_k4_scal)
      if (co_str_k4_scal /= str_k4_scal // 4_'          ') call print_and_register( 'get scalar kind=4 to kind=4 self failed.')

      str_k4_scal = co_str_k1_scal[1]
      print *, 4_'#' // str_k4_scal // 4_'#, len:', len(str_k4_scal)
      if (co_str_k4_scal /= str_k4_scal // 4_'          ') call print_and_register( 'get scalar kind=1 to kind=4 self failed.')

      str_k1_scal = co_str_k4_scal[1]
      print *, '#' // str_k1_scal // '#, len:', len(str_k1_scal)
      if (co_str_k1_scal /= str_k1_scal // '          ') call print_and_register( 'get scalar kind=4 to kind=1 self failed.')
    end if

    ! Do the same for arrays but on image 2
    if (me == 2) then
      str_k1_arr(:) = co_str_k1_arr(:)[2]
      print *, '#' // str_k1_arr(:) // '#, len:', len(str_k1_arr(1))
      if (any(str_k1_arr /= ['abc  ', 'EFG  ', 'klm  ', 'NOP  '])) &
        call print_and_register( 'get array kind=1 to kind=1 self failed.')

      print *, str_k4_arr
      str_k4_arr(:) = co_str_k4_arr(:)[2]
      print *, 4_'#' // str_k4_arr(:) // 4_'#, len:', len(str_k4_arr(1))
      if (any(str_k4_arr /= [4_'abc  ', 4_'EFG  ', 4_'klm  ', 4_'NOP  '])) &
        call print_and_register( 'get array kind=4 to kind=4 self failed.')

      str_k4_arr(:) = co_str_k1_arr(:)[2]
      print *, 4_'#' // str_k4_arr(:) // 4_'#, len:', len(str_k4_arr(1))
      if (any(str_k4_arr /= [ 4_'abc  ', 4_'EFG  ', 4_'klm  ', 4_'NOP  '])) &
        call print_and_register( 'get array kind=1 to kind=4 self failed.')

      str_k1_arr(:) = co_str_k4_arr(:)[2]
      print *, '#' // str_k1_arr(:) // '#, len:', len(str_k1_arr(1))
      if (any(str_k1_arr /= ['abc  ', 'EFG  ', 'klm  ', 'NOP  '])) &
        call print_and_register( 'get array kind=4 to kind=1 self failed.')
    end if

    sync all
    if (me == 1) then
      str_k1_scal = co_str_k1_scal[2]
      print *, '#' // str_k1_scal // '#, len:', len(str_k1_scal)
      if (co_str_k1_scal /= str_k1_scal // '          ') call print_and_register( 'get kind=1 to kind=1 image 2 failed.')

      str_k4_scal = co_str_k4_scal[2]
      print *, 4_'#' // str_k4_scal // 4_'#, len:', len(str_k4_scal)
      if (co_str_k4_scal /= str_k4_scal // 4_'          ') call print_and_register( 'get kind=4 to kind=4 image 2 failed.')
    else if (me == 2) then
      str_k4_scal = co_str_k1_scal[1]
      print *, 4_'#' // str_k4_scal // 4_'#, len:', len(str_k4_scal)
      if (co_str_k4_scal /= str_k4_scal // 4_'          ') call print_and_register( 'get kind=1 to kind=4 from image 1 failed.')

      str_k1_scal = co_str_k4_scal[1]
      print *, '#' // str_k1_scal // '#, len:', len(str_k1_scal)
      if (co_str_k1_scal /= str_k1_scal // '          ') call print_and_register( 'get kind=4 to kind=1 from image 1 failed.')
    end if

    str_k1_arr(:) = '#####'
    str_k4_arr(:) = 4_'#####'

    sync all

    if (me == 1) then
      str_k1_arr(1:2) = co_str_k1_arr(::2)[2]
      print *, str_k1_arr
      if (any(str_k1_arr /= ['abc  ', 'klm  ', '#####', '#####'])) &
        & call print_and_register( "strided get char arr kind 1 to kind 1 failed.")

      str_k4_arr(1:2) = co_str_k4_arr(::2)[2]
      print *, str_k4_arr
      if (any(str_k4_arr /= [4_'abc  ', 4_'klm  ', 4_'#####', 4_'#####'] )) &
        & call print_and_register( "strided get char arr kind 4 to kind 4 failed.")
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

end program get_convert_char_array

! vim:ts=2:sts=2:sw=2:

