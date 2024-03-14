program example

    use pm_kind, only: IK, LK, SK, RKD
    use pm_io, only: display_type
    use pm_timer, only: timer_type

    implicit none

    class(timer_type), allocatable :: timer
    logical(LK) :: failed
    integer :: i

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("timer = timer_type()")
                    timer = timer_type()
    call showtimerComponents()
    call disp%skip()

    call disp%skip()
    call disp%show("timer%clock = timer%time()")
                    timer%clock = timer%time()
    call showtimerComponents()
    call disp%skip()

    call disp%skip()
    call disp%show("timer%delta = timer%time(since = timer%start)")
                    timer%delta = timer%time(since = timer%start)
    call showtimerComponents()
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Idle for a certain number of seconds.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("timer = timer_type()")
                    timer = timer_type()
    call disp%show("call timer%wait(seconds = 0.1_RKD)")
                    call timer%wait(seconds = 0.1_RKD)
    call disp%show("timer%delta = timer%time() - timer%start")
                    timer%delta = timer%time() - timer%start
    call showtimerComponents()
    call disp%skip()

contains

    impure subroutine showtimerComponents()
        call disp%show("timer%start")
        call disp%show( timer%start )
        call disp%show("timer%clock")
        call disp%show( timer%clock )
        call disp%show("timer%delta")
        call disp%show( timer%delta )
        call disp%show("timer%resol")
        call disp%show( timer%resol )
    end subroutine

end program example