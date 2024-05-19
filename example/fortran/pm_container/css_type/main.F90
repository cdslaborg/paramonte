program example
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_arrayRange, only: getRange
    use pm_str, only: getCharVec

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block

        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type), allocatable :: scalar, vector(:), matrix(:,:), cube(:,:,:)

        call disp%skip()
        call disp%show("scalar = css_type(SKG_'ParaMonte')")
                        scalar = css_type(SKG_'ParaMonte')
        call disp%show("scalar")
        call disp%show( scalar , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("vector = css_type([character(3,SKG) :: 'a', 'aa', 'aaa'])")
                        vector = css_type([character(3,SKG) :: 'a', 'aa', 'aaa'])
        call disp%show("vector")
        call disp%show( vector , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("vector = css_type([character(3,SKG) :: 'a', 'aa', 'aaa'], trimmed = .false._LK)")
                        vector = css_type([character(3,SKG) :: 'a', 'aa', 'aaa'], trimmed = .false._LK)
        call disp%show("vector")
        call disp%show( vector , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("vector = css_type([character(3,SKG) :: 'a', 'aa', 'aaa'], trimmed = .true._LK) ! This will avoid the automatic application of `trim()` to the input string vector.")
                        vector = css_type([character(3,SKG) :: 'a', 'aa', 'aaa'], trimmed = .true._LK)
        call disp%show("vector")
        call disp%show( vector , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("matrix = css_type( reshape(getCharVec(getRange(SKG_'A', SKG_'X')), shape = [3, 8]) )")
                        matrix = css_type( reshape(getCharVec(getRange(SKG_'A', SKG_'X')), shape = [3, 8]) )
        call disp%show("matrix")
        call disp%show( matrix , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("cube = css_type( reshape(getCharVec(getRange(SKG_'A', SKG_'X')), shape = [2, 3, 4]) )")
                        cube = css_type( reshape(getCharVec(getRange(SKG_'A', SKG_'X')), shape = [2, 3, 4]) )
        call disp%show("cube")
        call disp%show( cube , deliml = SK_"""" )
        call disp%skip()

    end block

end program example