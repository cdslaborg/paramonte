program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: isValidDateTime, getYear

    implicit none

    integer :: i
    integer(IK) :: Values(9)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isValidDateTime(0_IK)")
    call disp%show( isValidDateTime(0_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getYear()")
    call disp%show( getYear() )
    call disp%show("isValidDateTime(getYear())")
    call disp%show( isValidDateTime(getYear()) )
    call disp%skip()

    call disp%skip()
    call disp%show("call date_and_time(values = Values(1:8))")
                    call date_and_time(values = Values(1:8))
    call disp%show("Values(1:8)")
    call disp%show( Values(1:8) )
    call disp%skip()
    do i = size(Values), 0, -1
        call disp%show("i")
        call disp%show( i )
        call disp%show("isValidDateTime(Values(1:i))")
        call disp%show( isValidDateTime(Values(1:i)) )
        call disp%skip()
    end do
    call disp%skip()

    call disp%skip()
    call disp%show("isValidDateTime(1999_IK, 2_IK, 28_IK) ! February is 28 days in non-leap years.")
    call disp%show( isValidDateTime(1999_IK, 2_IK, 28_IK) )
    call disp%show("isValidDateTime(1999_IK, 2_IK, 29_IK) ! February is 28 days in non-leap years.")
    call disp%show( isValidDateTime(1999_IK, 2_IK, 29_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isValidDateTime(2000_IK, 2_IK, 28_IK) ! February is 29 days in leap years.")
    call disp%show( isValidDateTime(2000_IK, 2_IK, 28_IK) )
    call disp%show("isValidDateTime(2000_IK, 2_IK, 29_IK) ! February is 29 days in leap years.")
    call disp%show( isValidDateTime(2000_IK, 2_IK, 29_IK) )
    call disp%skip()

end program example