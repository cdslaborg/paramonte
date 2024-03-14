program example

    use pm_kind, only: SK, LK, IK
    use pm_arrayMinMax, only: getMinMaxVal
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = SK_"main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("!character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => SK ! All kinds are supported.
        character(:,TKC), allocatable :: array
        character(1,TKC) :: minmax(2)
        call disp%skip()
        call disp%show("array = 'ParaMonte is a Machine Learning Library.'")
                        array = 'ParaMonte is a Machine Learning Library.'
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%show("minmax = getMinMaxVal(array)")
                        minmax = getMinMaxVal(array)
        call disp%show("minmax")
        call disp%show( minmax , deliml = SK_"""" )
        call disp%show("minmax = getMinMaxVal(array(1:0))")
                        minmax = getMinMaxVal(array(1:0))
        call disp%show("minmax")
        call disp%show( minmax , deliml = SK_"""" )
        call disp%show("minmax(1) == char(0)")
        call disp%show( minmax(1) == char(0) )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%")
    call disp%show("!character array.")
    call disp%show("!%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => SK ! All kinds are supported.
        character(10,TKC), allocatable :: array(:)
        character(10,TKC) :: minmax(2)
        call disp%skip()
        call disp%show('array = [character(10,TKC) :: "ParaMonte", "is", "a", "Monte", "Carlo", "Library."]')
                        array = [character(10,TKC) :: "ParaMonte", "is", "a", "Monte", "Carlo", "Library."]
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%show("minmax = getMinMaxVal(array)")
                        minmax = getMinMaxVal(array)
        call disp%show("minmax")
        call disp%show( minmax , deliml = SK_"""" )
        call disp%show("minmax = getMinMaxVal(array(1:0))")
                        minmax = getMinMaxVal(array(1:0))
        call disp%show("minmax")
        call disp%show( minmax , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("!string array.")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => SK ! All kinds are supported.
        use pm_container, only: css_type
        type(css_type), allocatable :: array(:)
        type(css_type) :: minmax(2)
        call disp%skip()
        call disp%show('array = css_type([character(10,TKC) :: "ParaMonte", "is", "a", "Monte", "Carlo", "Library."])')
                        array = css_type([character(10,TKC) :: "ParaMonte", "is", "a", "Monte", "Carlo", "Library."])
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%show("minmax = getMinMaxVal(array)")
                        minmax = getMinMaxVal(array)
        call disp%show("minmax")
        call disp%show( minmax , deliml = SK_"""" )
        call disp%show("minmax = getMinMaxVal(array(1:0))")
                        minmax = getMinMaxVal(array(1:0))
        call disp%show("[allocated(minmax(1)%val), allocated(minmax(2)%val)]")
        call disp%show( [allocated(minmax(1)%val), allocated(minmax(2)%val)] )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!integer array.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => IK ! All kinds are supported.
        integer(TKC), allocatable :: array(:)
        integer(TKC) :: minmax(2)
        call disp%skip()
        call disp%show('array = getUnifRand(-9, 9, getUnifRand(3_IK, 9_IK))')
                        array = getUnifRand(-9, 9, getUnifRand(3_IK, 9_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("minmax = getMinMaxVal(array)")
                        minmax = getMinMaxVal(array)
        call disp%show("minmax")
        call disp%show( minmax )
        call disp%show("minmax = getMinMaxVal(array(1:0))")
                        minmax = getMinMaxVal(array(1:0))
        call disp%show("minmax")
        call disp%show( minmax )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!logical array.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => LK ! All kinds are supported.
        logical(TKC), allocatable :: array(:)
        logical(TKC) :: minmax(2)
        call disp%skip()
        call disp%show('array = getUnifRand(.false., .true., getUnifRand(3_IK, 9_IK))')
                        array = getUnifRand(.false., .true., getUnifRand(3_IK, 9_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("minmax = getMinMaxVal(array)")
                        minmax = getMinMaxVal(array)
        call disp%show("minmax")
        call disp%show( minmax )
        call disp%show("minmax = getMinMaxVal(array(1:0))")
                        minmax = getMinMaxVal(array(1:0))
        call disp%show("minmax")
        call disp%show( minmax )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!complex array.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => CKS ! All kinds are supported.
        complex(TKC), allocatable :: array(:)
        complex(TKC) :: minmax(2)
        call disp%skip()
        call disp%show('array = getUnifRand((-9., -9.), (9., 9.), getUnifRand(3_IK, 9_IK))')
                        array = getUnifRand((-9., -9.), (9., 9.), getUnifRand(3_IK, 9_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("minmax = getMinMaxVal(array)")
                        minmax = getMinMaxVal(array)
        call disp%show("minmax")
        call disp%show( minmax )
        call disp%show("minmax = getMinMaxVal(array(1:0))")
                        minmax = getMinMaxVal(array(1:0))
        call disp%show("minmax")
        call disp%show( minmax )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!real array.")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! All kinds are supported.
        real(TKC), allocatable :: array(:)
        real(TKC) :: minmax(2)
        call disp%skip()
        call disp%show('array = getUnifRand(-9, 9, getUnifRand(3_IK, 9_IK))')
                        array = getUnifRand(-9, 9, getUnifRand(3_IK, 9_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("minmax = getMinMaxVal(array)")
                        minmax = getMinMaxVal(array)
        call disp%show("minmax")
        call disp%show( minmax )
        call disp%show("minmax = getMinMaxVal(array(1:0))")
                        minmax = getMinMaxVal(array(1:0))
        call disp%show("minmax")
        call disp%show( minmax )
        call disp%skip()
    end block

end program example