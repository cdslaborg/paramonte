program example

    use pm_kind, only: LK, SK
    use pm_logicalCompare, only: operator(>=)
    use pm_distUnif, only: setUnifRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! `.true.` evaluates to 1 and `.false.` evaluates to `0.`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show(".false._LK >= .false._LK")
    call disp%show( .false._LK >= .false._LK )
    call disp%skip()

    call disp%show(".false._LK >= .true._LK")
    call disp%show( .false._LK >= .true._LK )
    call disp%skip()

    call disp%show(".true._LK >= .false._LK")
    call disp%show( .true._LK >= .false._LK )
    call disp%skip()

    call disp%show(".true._LK >= .true._LK")
    call disp%show( .true._LK >= .true._LK )
    call disp%skip()

    call disp%show("[.true._LK, .false._LK] >= .true._LK")
    call disp%show( [.true._LK, .false._LK] >= .true._LK )
    call disp%skip()

    call disp%show("[.true._LK, .true._LK, .false._LK, .false._LK] >= [.true._LK, .false._LK, .true._LK, .false._LK]")
    call disp%show( [.true._LK, .true._LK, .false._LK, .false._LK] >= [.true._LK, .false._LK, .true._LK, .false._LK] )
    call disp%skip()

end program example