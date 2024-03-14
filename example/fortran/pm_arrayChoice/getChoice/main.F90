program example

    use pm_kind, only: LK ! All kinds are supported.
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayRange, only: getRange
    use pm_distUnif, only: getUnifRand
    use pm_arrayChoice, only: getChoice

    implicit none

    integer(IK) :: count, itry, ntry = 10
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

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
            call disp%show("choice = getChoice(array)")
                            choice = getChoice(array)
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )
            call disp%show("count = getUnifRand(1, len(array))")
                            count = getUnifRand(1, len(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.")
                            choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )
            call disp%show("choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.
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
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )
            call disp%show("count = getUnifRand(1, len(array))")
                            count = getUnifRand(1, len(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.")
                            choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice, deliml = SK_"""" )
            call disp%show("choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.
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
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, size(array))")
                            count = getUnifRand(1, size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.")
                            choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.
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
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, size(array))")
                            count = getUnifRand(1, size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.")
                            choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.
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
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, size(array))")
                            count = getUnifRand(1, size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.")
                            choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.
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
            call disp%show("choice = [getChoice(array)]")
                            choice = [getChoice(array)]
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("count = getUnifRand(1, size(array))")
                            count = getUnifRand(1, size(array))
            call disp%show("count")
            call disp%show( count )
            call disp%show("choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.")
                            choice = getChoice(array, count) ! draw randomly only `count` elements with replacement.
            call disp%show("choice")
            call disp%show( choice )
            call disp%show("choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.")
                            choice = getChoice(array, count, unique = .true._LK) ! draw randomly only `count` elements without replacement.
            call disp%show("choice")
            call disp%show( choice )

        end do
    end block

end program example