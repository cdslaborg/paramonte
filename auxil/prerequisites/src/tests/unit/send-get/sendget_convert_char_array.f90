!! Thoroughly test sendget, i.e. foo [N] = bar [M] in all variants
!!
!! Do simple tests for sendget(). These test comprise
!!
!! FOO = BAR [N]
!!
!! where
!!
!!  FOO                BAR                 images
!! character(len=20) character(len=10)   N == M == me
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
!! all of the above but for              N != M != me
!!
!! And may be some other, I've forgotten.
!!
!! Author: Andre Vehreschild, 2017


program sendget_convert_char_array

  implicit none

  character(kind=1, len=10), codimension[*] :: co_str_k1_src_scal
  character(kind=1, len=20), codimension[*] :: co_str_k1_dst_scal
  character(kind=4, len=10), codimension[*] :: co_str_k4_src_scal
  character(kind=4, len=20), codimension[*] :: co_str_k4_dst_scal

  character(kind=1, len=5), codimension[*] :: co_str_k1_src_arr(1:4)
  character(kind=1, len=7), codimension[*] :: co_str_k1_dst_arr(1:4)
  character(kind=4, len=5), codimension[*] :: co_str_k4_src_arr(1:4)
  character(kind=4, len=7), codimension[*] :: co_str_k4_dst_arr(1:4)

  logical :: error_printed=.false.

  associate(me => this_image(), np => num_images())
    if (np < 3) error stop 'Can not run with less than 3 images.'

    co_str_k1_src_scal = 'abcdefghij'
    co_str_k4_src_scal = 4_'abcdefghij'

    co_str_k1_src_arr(:) = ['abc', 'EFG', 'klm', 'NOP']
    co_str_k4_src_arr(:) = [4_'abc', 4_'EFG', 4_'klm', 4_'NOP']

    ! First check get/copy to self
    if (me == 1) then
      co_str_k1_dst_scal[1] = co_str_k1_src_scal[1]
      print *, '#' // co_str_k1_dst_scal // '#, len:', len(co_str_k1_dst_scal)
      if (co_str_k1_dst_scal /= co_str_k1_src_scal // '          ') error stop 'sendget scalar kind=1 to kind=1 self failed.'

      co_str_k4_dst_scal[1] = co_str_k4_src_scal[1]
      print *, 4_'#' // co_str_k4_dst_scal // 4_'#, len:', len(co_str_k4_dst_scal)
      if (co_str_k4_dst_scal /= co_str_k4_src_scal // 4_'          ') error stop 'sendget scalar kind=4 to kind=4 self failed.'

      co_str_k4_dst_scal[1] = co_str_k1_src_scal[1]
      print *, 4_'#' // co_str_k4_dst_scal // 4_'#, len:', len(co_str_k4_dst_scal)
      if (co_str_k4_dst_scal /= co_str_k4_src_scal // 4_'          ') error stop 'sendget scalar kind=1 to kind=4 self failed.'

      co_str_k1_dst_scal[1] = co_str_k4_src_scal[1]
      print *, '#' // co_str_k1_dst_scal // '#, len:', len(co_str_k1_dst_scal)
      if (co_str_k1_dst_scal /= co_str_k1_src_scal // '          ') error stop 'sendget scalar kind=4 to kind=1 self failed.'
    end if

    ! Do the same for arrays but on image 2
    if (me == 2) then
      co_str_k1_dst_arr(:)[2] = co_str_k1_src_arr(:)[2]
      print *, '#' // co_str_k1_dst_arr(:) // '#, len:', len(co_str_k1_dst_arr(1))
      if (any(co_str_k1_dst_arr /= ['abc    ', 'EFG    ', 'klm    ', 'NOP    '])) &
        & error stop 'sendget array kind=1 to kind=1 self failed.'

      co_str_k4_dst_arr(:)[2] = co_str_k4_src_arr(:)[2]
      print *, 4_'#' // co_str_k4_dst_arr(:) // 4_'#, len:', len(co_str_k4_dst_arr(1))
      if (any(co_str_k4_dst_arr /= [4_'abc    ', 4_'EFG    ', 4_'klm    ', 4_'NOP    '])) &
        & error stop 'sendget array kind=4 to kind=4 self failed.'

      co_str_k4_dst_arr(:)[2] = co_str_k1_src_arr(:)[2]
      print *, 4_'#' // co_str_k4_dst_arr(:) // 4_'#, len:', len(co_str_k4_dst_arr(1))
      if (any(co_str_k4_dst_arr /= [ 4_'abc    ', 4_'EFG    ', 4_'klm    ', 4_'NOP    '])) &
        & error stop 'sendget array kind=1 to kind=4 self failed.'

      co_str_k1_dst_arr(:)[2] = co_str_k4_src_arr(:)[2]
      print *, '#' // co_str_k1_dst_arr(:) // '#, len:', len(co_str_k1_dst_arr(1))
      if (any(co_str_k1_dst_arr /= ['abc    ', 'EFG    ', 'klm    ', 'NOP    '])) &
        & error stop 'sendget array kind=4 to kind=1 self failed.'
    end if

    co_str_k1_dst_arr(:) = '#######'
    co_str_k4_dst_arr(:) = 4_'#######'

    sync all
    if (me == 2) then
      co_str_k1_dst_scal[3] = co_str_k1_src_scal[1]

      co_str_k4_dst_scal[3] = co_str_k4_src_scal[1]

      co_str_k1_dst_arr(:)[3] = co_str_k1_src_arr(:)[1]

      co_str_k4_dst_arr(:)[3] = co_str_k4_src_arr(:)[1]
    end if

    sync all
    if (me == 3) then
      print *, '#' // co_str_k1_dst_scal // '#, len:', len(co_str_k1_dst_scal)
      if (co_str_k1_dst_scal /= co_str_k1_src_scal // '          ') &
        & error stop 'sendget kind=1 to kind=1 from image 1 to image 3 failed.'

      print *, 4_'#' // co_str_k4_dst_scal // 4_'#, len:', len(co_str_k4_dst_scal)
      if (co_str_k4_dst_scal /= co_str_k4_src_scal // 4_'          ') &
        & error stop 'sendget kind=4 to kind=4 from image 1 to image 3 failed.'

      print *, '#' // co_str_k1_dst_arr // '#, len:', len(co_str_k1_dst_arr(1))
      if (any(co_str_k1_dst_arr /= [co_str_k1_src_arr(1) // '  ', co_str_k1_src_arr(2) // '  ', &
        &                           co_str_k1_src_arr(3) // '  ', co_str_k1_src_arr(4) // '  '])) &
        & error stop 'sendget array kind=1 to kind=1 from image 1 to image 3 failed.'

      print *, 4_'#' // co_str_k4_dst_arr // 4_'#, len:', len(co_str_k4_dst_arr(1))
      if (any(co_str_k4_dst_arr /= [co_str_k4_src_arr(1) // 4_'  ', co_str_k4_src_arr(2) // 4_'  ', &
        &                           co_str_k4_src_arr(3) // 4_'  ', co_str_k4_src_arr(4) // 4_'  '])) &
        & error stop 'sendget array kind=4 to kind=4 from image 1 to image 3 failed.'
    end if

    co_str_k1_dst_arr(:) = '#######'
    co_str_k4_dst_arr(:) = 4_'#######'

    sync all
    if (me == 2) then
      co_str_k1_dst_scal[3] = co_str_k4_src_scal[1]

      co_str_k4_dst_scal[3] = co_str_k1_src_scal[1]

      co_str_k1_dst_arr(:)[3] = co_str_k4_src_arr(:)[1]

      co_str_k4_dst_arr(:)[3] = co_str_k1_src_arr(:)[1]
    end if

    sync all
    if (me == 3) then
      print *, '#' // co_str_k1_dst_scal // '#, len:', len(co_str_k1_dst_scal)
      if (co_str_k1_dst_scal /= co_str_k1_src_scal // '          ') &
        & error stop 'sendget kind=4 to kind=1 from image 1 to image 3 failed.'

      print *, 4_'#' // co_str_k4_dst_scal // 4_'#, len:', len(co_str_k4_dst_scal)
      if (co_str_k4_dst_scal /= co_str_k4_src_scal // 4_'          ') &
        & error stop 'sendget kind=1 to kind=4 from image 1 to image 3 failed.'

      print *, '#' // co_str_k1_dst_arr // '#, len:', len(co_str_k1_dst_arr(1))
      if (any(co_str_k1_dst_arr /= [co_str_k1_src_arr(1) // '  ', co_str_k1_src_arr(2) // '  ', &
        &                           co_str_k1_src_arr(3) // '  ', co_str_k1_src_arr(4) // '  '])) &
        & error stop 'sendget array kind=4 to kind=1 from image 1 to image 3 failed.'

      print *, 4_'#' // co_str_k4_dst_arr // 4_'#, len:', len(co_str_k4_dst_arr(1))
      if (any(co_str_k4_dst_arr /= [co_str_k4_src_arr(1) // 4_'  ', co_str_k4_src_arr(2) // 4_'  ', &
        &                           co_str_k4_src_arr(3) // 4_'  ', co_str_k4_src_arr(4) // 4_'  '])) &
        & error stop 'sendget array kind=1 to kind=4 from image 1 to image 3 failed.'
    end if

    sync all

    co_str_k1_dst_arr(:) = '#######'
    co_str_k4_dst_arr(:) = 4_'#######'

    ! Now strided.
    sync all
    if (me == 2) then
      co_str_k1_dst_arr(4:1:-2)[3] = co_str_k4_src_arr(::2)[1]

      co_str_k4_dst_arr(4:1:-2)[3] = co_str_k1_src_arr(::2)[1]
    end if

    sync all
    if (me == 3) then
      print *, '#' // co_str_k1_dst_arr // '#, len:', len(co_str_k1_dst_arr(1))
      if (any(co_str_k1_dst_arr /= ['#######', co_str_k1_src_arr(3) // '  ', &
        &                           '#######', co_str_k1_src_arr(1) // '  '])) &
        & error stop 'sendget strided kind=4 to kind=1 from image 1 to image 3 failed.'

      print *, 4_'#' // co_str_k4_dst_arr // 4_'#, len:', len(co_str_k4_dst_arr(1))
      if (any(co_str_k4_dst_arr /= [4_'#######', co_str_k4_src_arr(3) // 4_'  ', &
        &                           4_'#######', co_str_k4_src_arr(1) // 4_'  '])) &
        & error stop 'sendget strided kind=1 to kind=4 from image 1 to image 3 failed.'
    end if

    select case(me)
      case(1)
        if (error_printed) error stop
        sync images([2,3])
        print *, 'Test passed.'
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

end program sendget_convert_char_array

! vim:ts=2:sts=2:sw=2:
