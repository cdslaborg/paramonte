program example

    use pm_kind, only: SK, LK
    use pm_kind, only: IK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_arrayComplement, only: getComplementRange
    use pm_arraySort, only: setSorted
    use pm_arrayUnique, only: getUnique

    implicit none

    integer(IK)     , allocatable   :: SetA_IK(:)
    integer(IK)                     :: start, stop, step

    type(display_type)              :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the complement of integer vector A in B: B \ A.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    SetA_IK = [integer(IK) :: 4, 5, 4, 1, 2, 1, 5]

    call disp%skip()
    call disp%show("SetA_IK")
    call disp%show( SetA_IK )

    call disp%show("start = 0; stop = 3; step = 1")
                    start = 0; stop = 3; step = 1
    call disp%show("getComplementRange(SetA_IK, start, stop, step)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step) )

    call disp%show("start = 0; stop = 6; step = 2")
                    start = 0; stop = 6; step = 2
    call disp%show("getComplementRange(SetA_IK, start, stop, step)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step) )

    call disp%show("start = 0; stop = 6; step = 1")
                    start = 0; stop = 6; step = 1
    call disp%show("getComplementRange(SetA_IK, start, stop, step)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step) )

    call disp%show("start = 10; stop = 6; step = -1")
                    start = 10; stop = 6; step = -1
    call disp%show("getComplementRange(SetA_IK, start, stop, step)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setSorted(SetA_IK)")
                    call setSorted(SetA_IK)
    call disp%show("SetA_IK")
    call disp%show( SetA_IK )

    call disp%show("start = 0; stop = 3; step = 1")
                    start = 0; stop = 3; step = 1
    call disp%show("getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .false._LK)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .false._LK) )

    call disp%show("start = 0; stop = 6; step = 2")
                    start = 0; stop = 6; step = 2
    call disp%show("getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .false._LK)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .false._LK) )
    call disp%skip()

    call disp%skip()
    call disp%show("SetA_IK = getUnique(SetA_IK)")
                    SetA_IK = getUnique(SetA_IK)
    call disp%show("SetA_IK")
    call disp%show( SetA_IK )

    call disp%show("start = 0; stop = 3; step = 1")
                    start = 0; stop = 3; step = 1
    call disp%show("getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .true._LK)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .true._LK) )

    call disp%show("start = 0; stop = 6; step = 2")
                    start = 0; stop = 6; step = 2
    call disp%show("getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .true._LK)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .true._LK) )

    call disp%show("start = 10; stop = 15; step = 1")
                    start = 10; stop = 15; step = 1
    call disp%show("getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .true._LK)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .true._LK) )

    call disp%skip()

    SetA_IK = [integer(IK) :: 5, 4, 3, 1, 1]

    call disp%skip()
    call disp%show("call setSorted(SetA_IK)")
                    call setSorted(SetA_IK)
    call disp%show("SetA_IK")
    call disp%show( SetA_IK )

    call disp%show("start = 6; stop = -3; step = -2")
                    start = 6; stop = -1; step = -2
    call disp%show("getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .false._LK)")
    call disp%show( getComplementRange(SetA_IK, start, stop, step, sorted = .true._LK, unique = .false._LK) )
    call disp%skip()

end program example