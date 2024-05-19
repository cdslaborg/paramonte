program example

    use pm_kind, only: LK ! All kinds are supported.
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayRange, only: getRange
    use pm_arrayChoice, only: getChoice
    use pm_distUnif, only: getUnifRand
    use pm_arrayChange, only: getChange

    implicit none

    integer(IK) :: count, itry, ntry = 10
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        integer(IK) :: step
        character(:), allocatable :: start, stop, change

        do itry = 1, 2
            call disp%skip
            call disp%show("count = getChoice([1, 2]); start = 'A'; stop = 'A'; step = getChoice([-1, 1])")
                            count = getChoice([1, 2]); start = 'A'; stop = 'A'; step = getChoice([-1, 1])
            call disp%show("count")
            call disp%show( count )
            call disp%show("change = getChange(count, start, stop, step)")
                            change = getChange(count, start, stop, step)
            call disp%show("change")
            call disp%show( change, deliml = SK_"""" )
        end do

        do itry = 1, ntry

            call disp%skip
            call disp%show("count = getUnifRand(1, 10); start = 'A'; stop = 'F'; step = 1")
                            count = getUnifRand(1, 10); start = 'A'; stop = 'F'; step = 1
            call disp%show("count")
            call disp%show( count )
            call disp%show("change = getChange(count, start, stop, step)")
                            change = getChange(count, start, stop, step)
            call disp%show("change")
            call disp%show( change, deliml = SK_"""" )

            call disp%skip
            call disp%show("count = getUnifRand(1, 6); start = 'F'; stop = 'A'; step = -1")
                            count = getUnifRand(1, 6); start = 'F'; stop = 'A'; step = -1
            call disp%show("count")
            call disp%show( count )
            call disp%show("change = getChange(count, start, stop, step, unique = .true._LK)")
                            change = getChange(count, start, stop, step, unique = .true._LK)
            call disp%show("change")
            call disp%show( change, deliml = SK_"""" )

            call disp%skip
            call disp%show("count = getUnifRand(5, 10); start = 'A'; stop = 'z'; step = 3")
                            count = getUnifRand(5, 10); start = 'A'; stop = 'z'; step = 3
            call disp%show("count")
            call disp%show( count )
            call disp%show("change = getChange(count, start, stop, step, unique = .true._LK)")
                            change = getChange(count, start, stop, step, unique = .true._LK)
            call disp%show("change")
            call disp%show( change, deliml = SK_"""" )

        end do
    end block

    block
        use pm_kind, only: IKG => IK ! all integer kinds are supported.
        integer(IKG) :: step
        integer(IKG), allocatable :: start, stop, change(:)
        do itry = 1, ntry

            call disp%skip
            call disp%show("count = getUnifRand(1, 10); start = 0; stop = 9; step = 1")
                            count = getUnifRand(1, 10); start = 0; stop = 9; step = 1
            call disp%show("count")
            call disp%show( count )
            call disp%show("change = getChange(count, start, stop, step)")
                            change = getChange(count, start, stop, step)
            call disp%show("change")
            call disp%show( change, deliml = SK_"""" )

            call disp%skip
            call disp%show("count = getUnifRand(5, 10); start = 9; stop = -9; step = -2")
                            count = getUnifRand(5, 10); start = 9; stop = -9; step = -2
            call disp%show("count")
            call disp%show( count )
            call disp%show("change = getChange(count, start, stop, step, unique = .true._LK)")
                            change = getChange(count, start, stop, step, unique = .true._LK)
            call disp%show("change")
            call disp%show( change, deliml = SK_"""" )

            call disp%skip
            call disp%show("count = getUnifRand(5, 15); start = 5; stop = 1; step = -2")
                            count = getUnifRand(5, 15); start = 5; stop = 1; step = -2
            call disp%show("count")
            call disp%show( count )
            call disp%show("change = getChange(count, start, stop, step, unique = .false._LK)")
                            change = getChange(count, start, stop, step, unique = .false._LK)
            call disp%show("change")
            call disp%show( change, deliml = SK_"""" )

        end do
    end block

end program example