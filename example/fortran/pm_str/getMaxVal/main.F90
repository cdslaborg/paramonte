program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_arrayFind, only: getLoc
    use pm_str, only: getMaxVal

    implicit none

    logical(LK)     , allocatable :: Mask(:)
    character(:, SK), allocatable :: str
    integer :: i

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getMaxVal(SK_'ParaMonte')")
    call disp%show( getMaxVal(SK_'ParaMonte') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getMaxVal(SK_'The ParaMonte Parallel Library')")
    call disp%show( getMaxVal(SK_'The ParaMonte Parallel Library') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("str = SK_'The ParaMonte Parallel Library'")
                    str = SK_'The ParaMonte Parallel Library'
    call disp%show("Mask = [logical(LK) :: (.true., i = 1, len(str))]; Mask(getLoc(str, SK_'y')) = .false._LK ! set all `y` locations in the Mask to `.false.`")
                    Mask = [logical(LK) :: (.true., i = 1, len(str))]; Mask(getLoc(str, SK_'y')) = .false._LK ! set all `y` locations in the Mask to `.false.`
    call disp%show("Mask")
    call disp%show( Mask )
    call disp%show("getMaxVal(str, Mask)")
    call disp%show( getMaxVal(str, Mask) , deliml = SK_"""" )
    call disp%show("getMaxVal(str)")
    call disp%show( getMaxVal(str) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getMaxVal([character(11,SK) :: 'The Fortran', 'Programming', 'Language']) ! Note that 'Language' is padded with blanks.")
    call disp%show( getMaxVal([character(11,SK) :: 'The Fortran', 'Programming', 'Language']) , deliml = SK_"""" )
    call disp%skip()

end program example