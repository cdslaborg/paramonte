program example
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: display_type
    use pm_arrayRange, only: getRange
    use pm_str, only: getCharVec

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

#if PDT_ENABLED
    block

        use pm_kind, only: SKC => SK
        use pm_container, only: css_pdt
        type(css_pdt), allocatable :: scalar, vector(:), matrix(:,:), cube(:,:,:)

        call disp%skip()
        call disp%show("scalar = css_pdt(SKC_'ParaMonte')")
                        scalar = css_pdt(SKC_'ParaMonte')
        call disp%show("scalar")
        call disp%show( scalar , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("vector = css_pdt([character(3,SKC) :: 'a', 'aa', 'aaa'])")
                        vector = css_pdt([character(3,SKC) :: 'a', 'aa', 'aaa'])
        call disp%show("vector")
        call disp%show( vector , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("vector = css_pdt([character(3,SKC) :: 'a', 'aa', 'aaa'], trimmed = .false._LK)")
                        vector = css_pdt([character(3,SKC) :: 'a', 'aa', 'aaa'], trimmed = .false._LK)
        call disp%show("vector")
        call disp%show( vector , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("vector = css_pdt([character(3,SKC) :: 'a', 'aa', 'aaa'], trimmed = .true._LK) ! This will avoid the automatic application of `trim()` to the input string vector.")
                        vector = css_pdt([character(3,SKC) :: 'a', 'aa', 'aaa'], trimmed = .true._LK)
        call disp%show("vector")
        call disp%show( vector , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("matrix = css_pdt( reshape(getCharVec(getRange(SKC_'A', SKC_'X')), shape = [3, 8]) )")
                        matrix = css_pdt( reshape(getCharVec(getRange(SKC_'A', SKC_'X')), shape = [3, 8]) )
        call disp%show("matrix")
        call disp%show( matrix , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("cube = css_pdt( reshape(getCharVec(getRange(SKC_'A', SKC_'X')), shape = [2, 3, 4]) )")
                        cube = css_pdt( reshape(getCharVec(getRange(SKC_'A', SKC_'X')), shape = [2, 3, 4]) )
        call disp%show("cube")
        call disp%show( cube , deliml = SK_"""" )
        call disp%skip()

    end block
#else
    call disp%show("PDT examples are disabled in the current ParaMonte library configuration.")
#endif

end program example