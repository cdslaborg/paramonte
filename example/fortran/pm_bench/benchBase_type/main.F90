program example

    use pm_kind, only: IK, LK, SK, RKD
    use pm_io, only: display_type
    use pm_bench, only: benchBase_type
    use pm_timer, only: timerCPU_type
    use pm_timer, only: timerDAT_type
    use pm_timer, only: timerSYS_type

    implicit none

    type(benchBase_type) :: benchBase
    real(RKD)           :: unifrnd
    real(RKD)           :: sumTime = 0._RKD
    real(RKD)           :: sumUnif = 0._RKD
    integer             :: i

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("benchBase = benchBase_type(name = SK_'random_number')")
                    benchBase = benchBase_type(name = SK_'random_number')
    call showTimerComponents()
    call disp%skip()

    call disp%skip()
    call disp%show("allocate(benchBase%timing%values(1000))")
    call disp%show("do i = 1, 1000")
    call disp%show("    benchBase%timer%start = benchBase%timer%time()")
    call disp%show("    call random_number(unifrnd)")
    call disp%show("    benchBase%timing%values(i) = benchBase%timer%time(since = benchBase%timer%start)")
    call disp%show("    sumTime = sumTime + benchBase%timing%values(i)")
    call disp%show("    if (sumTime > benchBase%minsec) exit")
    call disp%show("    sumUnif = sumUnif + unifrnd")
    call disp%show("end do")
                    allocate(benchBase%timing%values(1000))
                    do i = 1, 1000
                        benchBase%timer%start = benchBase%timer%time()
                        call random_number(unifrnd)
                        benchBase%timing%values(i) = benchBase%timer%time(since = benchBase%timer%start)
                        sumTime = sumTime + benchBase%timing%values(i)
                        if (sumTime > benchBase%minsec) exit
                        sumUnif = sumUnif + unifrnd
                    end do
    call disp%show("sum(benchBase%timing%values) / size(benchBase%timing%values)")
    call disp%show( sum(benchBase%timing%values) / size(benchBase%timing%values) )
    call disp%skip()

contains

    impure subroutine showTimerComponents()
        call disp%show("benchBase%name")
        call disp%show( benchBase%name , deliml = SK_"""" )
        call disp%show("benchBase%minsec")
        call disp%show( benchBase%minsec )
        call disp%show("benchBase%miniter")
        call disp%show( benchBase%miniter )
    end subroutine

end program example