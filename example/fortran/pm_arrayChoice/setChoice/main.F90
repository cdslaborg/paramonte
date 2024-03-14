program example

    use pm_kind, only: LK ! All kinds are supported.
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayRange, only: getRange
    use pm_distUnif, only: getUnifRand
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: xoshiro256ssw_type
    use pm_arrayChoice, only: setChoice, getChoice

    implicit none

    type(display_type) :: disp
    integer(IK) :: count, itry, ntry = 10
    type(xoshiro256ssw_type) :: rng
    disp = display_type(file = "main.out.F90")
    rng = xoshiro256ssw_type()

    block
        character(:), allocatable :: choice, array
        do itry = 1, ntry

            call disp%skip
            call disp%show("count = getUnifRand(4, 10)")
                            count = getUnifRand(4, 10)
            call disp%show("array = getUnifRand(repeat('A', count), repeat('Z', count)) ! generate random array for illustration.")
                            array = getUnifRand(repeat('A', count), repeat('Z', count)) ! generate random array for illustration.
            call disp%show("array")
            call disp%show( array, deliml = SK_"""" )
            call disp%show("call setResized(choice, 1)")
                            call setResized(choice, 1)
            call disp%show("call setChoice(rng, choice, array)")
                            call setChoice(rng, choice, array)
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )
            call disp%show("count = getUnifRand(1, 2 * len(array))")
                            count = getUnifRand(1, 2 * len(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.")
                            call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )
            call disp%show("count = getUnifRand(1, len(array))")
                            count = getUnifRand(1, len(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )

        end do
    end block

    block
        character(2), allocatable :: choice(:), array(:)
        do itry = 1, ntry

            call disp%skip
            call disp%show("count = getUnifRand(4, 10)")
                            count = getUnifRand(4, 10)
            call disp%show("array = getUnifRand('AA', 'ZZ', count) ! generate random array for illustration.")
                            array = getUnifRand('AA', 'ZZ', count) ! generate random array for illustration.
            call disp%show("array")
            call disp%show( array, deliml = SK_"""" )
            call disp%show("call setResized(choice, 1)")
                            call setResized(choice, 1)
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )
            call disp%show("count = getUnifRand(1, len(array))")
                            count = getUnifRand(1, len(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.")
                            call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )
            call disp%show("count = getUnifRand(1, size(array))")
                            count = getUnifRand(1, size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )

        end do
    end block

    block
        integer, allocatable :: choice(:), array(:)
        do itry = 1, ntry

            call disp%skip
            call disp%show("array = getUnifRand(1, 20, getUnifRand(4, 10)) ! generate random array for illustration.")
                            array = getUnifRand(1, 20, getUnifRand(4, 10)) ! generate random array for illustration.
            call disp%show("array")
            call disp%show( array )
            call disp%show("call setResized(choice, 1)")
                            call setResized(choice, 1)
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, 2 * size(array))")
                            count = getUnifRand(1, 2 * size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.")
                            call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, size(array))")
                            count = getUnifRand(1, size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.
            call disp%show("choice")
            call disp%show( choice )

        end do
    end block

    block
        logical, allocatable :: choice(:), array(:)
        do itry = 1, ntry

            call disp%skip
            call disp%show("array = getUnifRand(.false., .true., getUnifRand(4, 10)) ! generate random array for illustration.")
                            array = getUnifRand(.false., .true., getUnifRand(4, 10)) ! generate random array for illustration.
            call disp%show("array")
            call disp%show( array )
            call disp%show("call setResized(choice, 1)")
                            call setResized(choice, 1)
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, 2 * size(array))")
                            count = getUnifRand(1, 2 * size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.")
                            call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, size(array))")
                            count = getUnifRand(1, size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.
            call disp%show("choice")
            call disp%show( choice )

        end do
    end block

    block
        complex, allocatable :: choice(:), array(:)
        do itry = 1, ntry

            call disp%skip
            call disp%show("array = getUnifRand((0., 0.), (1., 1.), getUnifRand(4, 10)) ! generate random array for illustration.")
                            array = getUnifRand((0., 0.), (1., 1.), getUnifRand(4, 10)) ! generate random array for illustration.
            call disp%show("array")
            call disp%show( array )
            call disp%show("call setResized(choice, 1)")
                            call setResized(choice, 1)
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, 2 * size(array))")
                            count = getUnifRand(1, 2 * size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.")
                            call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, size(array))")
                            count = getUnifRand(1, size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.
            call disp%show("choice")
            call disp%show( choice )

        end do
    end block

    block
        real, allocatable :: choice(:), array(:)
        do itry = 1, ntry

            call disp%skip
            call disp%show("array = getUnifRand(0., 1., getUnifRand(4, 10)) ! generate random array for illustration.")
                            array = getUnifRand(0., 1., getUnifRand(4, 10)) ! generate random array for illustration.
            call disp%show("array")
            call disp%show( array )
            call disp%show("call setResized(choice, 1)")
                            call setResized(choice, 1)
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, 2 * size(array))")
                            count = getUnifRand(1, 2 * size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.")
                            call setChoice(rng, choice, array) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, size(array))")
                            count = getUnifRand(1, size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("call setResized(choice, count)")
                            call setResized(choice, count)
            call disp%show("call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            call setChoice(rng, choice, array, unique = .true._LK) ! draw randomly only `count` elements without replacement.
            call disp%show("choice")
            call disp%show( choice )

        end do
    end block

end program example