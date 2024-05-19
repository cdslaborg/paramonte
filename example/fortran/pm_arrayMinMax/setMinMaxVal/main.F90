program example

    use pm_kind, only: SK, LK, IK
    use pm_arrayMinMax, only: setMinMaxVal
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
        use pm_kind, only: TKG => SK ! All kinds are supported.
        character(:,TKG), allocatable :: array
        character(1,TKG) :: vmin, vmax
        call disp%skip()
        call disp%show("array = 'ParaMonte is a Machine Learning Library.'")
                        array = 'ParaMonte is a Machine Learning Library.'
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%show("call setMinMaxVal(array, vmin, vmax)")
                        call setMinMaxVal(array, vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] , deliml = SK_"""" )
        call disp%show("call setMinMaxVal(array(1:0), vmin, vmax)")
                        call setMinMaxVal(array(1:0), vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] , deliml = SK_"""" )
        call disp%show("vmin == char(0)")
        call disp%show( vmin == char(0) )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%")
    call disp%show("!character array.")
    call disp%show("!%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => SK ! All kinds are supported.
        character(10,TKG), allocatable :: array(:)
        character(10,TKG) :: vmin, vmax
        call disp%skip()
        call disp%show('array = [character(10,TKG) :: "ParaMonte", "is", "a", "Monte", "Carlo", "Library."]')
                        array = [character(10,TKG) :: "ParaMonte", "is", "a", "Monte", "Carlo", "Library."]
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%show("call setMinMaxVal(array, vmin, vmax)")
                        call setMinMaxVal(array, vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] , deliml = SK_"""" )
        call disp%show("call setMinMaxVal(array(1:0), vmin, vmax)")
                        call setMinMaxVal(array(1:0), vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("!string array.")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => SK ! All kinds are supported.
        use pm_container, only: css_type
        type(css_type), allocatable :: array(:)
        type(css_type) :: vmin, vmax
        call disp%skip()
        call disp%show('array = css_type([character(10,TKG) :: "ParaMonte", "is", "a", "Monte", "Carlo", "Library."])')
                        array = css_type([character(10,TKG) :: "ParaMonte", "is", "a", "Monte", "Carlo", "Library."])
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%show("call setMinMaxVal(array, vmin, vmax)")
                        call setMinMaxVal(array, vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] , deliml = SK_"""" )
        call disp%show("call setMinMaxVal(array(1:0), vmin, vmax)")
                        call setMinMaxVal(array(1:0), vmin, vmax)
        call disp%show("[allocated(vmin%val), allocated(vmax%val)]")
        call disp%show( [allocated(vmin%val), allocated(vmax%val)] )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!integer array.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => IK ! All kinds are supported.
        integer(TKG), allocatable :: array(:)
        integer(TKG) :: vmin, vmax
        call disp%skip()
        call disp%show('array = getUnifRand(-9, 9, getUnifRand(3_IK, 9_IK))')
                        array = getUnifRand(-9, 9, getUnifRand(3_IK, 9_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("call setMinMaxVal(array, vmin, vmax)")
                        call setMinMaxVal(array, vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] )
        call disp%show("call setMinMaxVal(array(1:0), vmin, vmax)")
                        call setMinMaxVal(array(1:0), vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!logical array.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => LK ! All kinds are supported.
        logical(TKG), allocatable :: array(:)
        logical(TKG) :: vmin, vmax
        call disp%skip()
        call disp%show('array = getUnifRand(.false., .true., getUnifRand(3_IK, 9_IK))')
                        array = getUnifRand(.false., .true., getUnifRand(3_IK, 9_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("call setMinMaxVal(array, vmin, vmax)")
                        call setMinMaxVal(array, vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] )
        call disp%show("call setMinMaxVal(array(1:0), vmin, vmax)")
                        call setMinMaxVal(array(1:0), vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!complex array.")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => CKS ! All kinds are supported.
        complex(TKG), allocatable :: array(:)
        complex(TKG) :: vmin, vmax
        call disp%skip()
        call disp%show('array = getUnifRand((-9., -9.), (9., 9.), getUnifRand(3_IK, 9_IK))')
                        array = getUnifRand((-9., -9.), (9., 9.), getUnifRand(3_IK, 9_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("call setMinMaxVal(array, vmin, vmax)")
                        call setMinMaxVal(array, vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] )
        call disp%show("call setMinMaxVal(array(1:0), vmin, vmax)")
                        call setMinMaxVal(array(1:0), vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!real array.")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        real(TKG), allocatable :: array(:)
        real(TKG) :: vmin, vmax
        call disp%skip()
        call disp%show('array = getUnifRand(-9, 9, getUnifRand(3_IK, 9_IK))')
                        array = getUnifRand(-9, 9, getUnifRand(3_IK, 9_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("call setMinMaxVal(array, vmin, vmax)")
                        call setMinMaxVal(array, vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] )
        call disp%show("call setMinMaxVal(array(1:0), vmin, vmax)")
                        call setMinMaxVal(array(1:0), vmin, vmax)
        call disp%show("[vmin, vmax]")
        call disp%show( [vmin, vmax] )
        call disp%skip()
    end block

end program example