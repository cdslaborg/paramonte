program example

    use pm_kind, only: SK
    use pm_kind, only: CK ! All other kinds are also supported.
    use pm_complexCompareAll, only: operator(<=)
    use pm_distUnif, only: setUnifRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Verify if all `complex` components follow the relation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("(+0._CK, +0._CK) <= (+0._CK, +0._CK)")
    call disp%show( (+0._CK, +0._CK) <= (+0._CK, +0._CK) )
    call disp%skip()

    call disp%show("(+0._CK, +0._CK) <= (+1._CK, +0._CK)")
    call disp%show( (+0._CK, +0._CK) <= (+1._CK, +0._CK) )
    call disp%skip()

    call disp%show("(+0._CK, +0._CK) <= (+0._CK, +1._CK)")
    call disp%show( (+0._CK, +0._CK) <= (+0._CK, +1._CK) )
    call disp%skip()

    call disp%show("(+0._CK, +0._CK) <= (-1._CK, +0._CK)")
    call disp%show( (+0._CK, +0._CK) <= (-1._CK, +0._CK) )
    call disp%skip()

    call disp%show("(+0._CK, +0._CK) <= (+0._CK, -1._CK)")
    call disp%show( (+0._CK, +0._CK) <= (+0._CK, -1._CK) )
    call disp%skip()

    call disp%show("(+0._CK, +0._CK) <= (+1._CK, -1._CK)")
    call disp%show( (+0._CK, +0._CK) <= (+1._CK, -1._CK) )
    call disp%skip()

    call disp%show("(+0._CK, +0._CK) <= (-1._CK, +1._CK)")
    call disp%show( (+0._CK, +0._CK) <= (-1._CK, +1._CK) )
    call disp%skip()

    call disp%show("(+0._CK, +0._CK) <= (+1._CK, +1._CK)")
    call disp%show( (+0._CK, +0._CK) <= (+1._CK, +1._CK) )
    call disp%skip()

    call disp%show("(+0._CK, +0._CK) <= [(+0._CK, +0._CK), (+1._CK, +0._CK), (+0._CK, +1._CK), (-1._CK, +0._CK), (+1._CK, -1._CK), (-1._CK, +1._CK), (-1._CK, -1._CK), (+1._CK, +1._CK)]")
    call disp%show( (+0._CK, +0._CK) <= [(+0._CK, +0._CK), (+1._CK, +0._CK), (+0._CK, +1._CK), (-1._CK, +0._CK), (+1._CK, -1._CK), (-1._CK, +1._CK), (-1._CK, -1._CK), (+1._CK, +1._CK)] )
    call disp%skip()

end program example