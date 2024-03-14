program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_arrayFind, only: getLoc
    use pm_str, only: getMinVal

    implicit none

    logical(LK)     , allocatable :: Mask(:)
    character(:, SK), allocatable :: str
    integer :: i

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getMinVal(SK_'ParaMonte')")
    call disp%show( getMinVal(SK_'ParaMonte') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getMinVal(SK_'The ParaMonte Parallel Library')")
    call disp%show( getMinVal(SK_'The ParaMonte Parallel Library') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("str = SK_'The ParaMonte Parallel Library'")
                    str = SK_'The ParaMonte Parallel Library'
    call disp%show("Mask = [logical(LK) :: (.true., i = 1, len(str))]; Mask(getLoc(str, SK_' ')) = .false._LK ! set all whitespace locations in the Mask to `.false.`")
                    Mask = [logical(LK) :: (.true., i = 1, len(str))]; Mask(getLoc(str, SK_' ')) = .false._LK ! set all whitespace locations in the Mask to `.false.`
    call disp%show("Mask")
    call disp%show( Mask )
    call disp%show("getMinVal(str, Mask)")
    call disp%show( getMinVal(str, Mask) , deliml = SK_"""" )
    call disp%show("getMinVal(str)")
    call disp%show( getMinVal(str) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getMinVal([character(11,SK) :: 'The Fortran', 'Programming', 'Language']) ! Note that 'Language' is padded with blanks.")
    call disp%show( getMinVal([character(11,SK) :: 'The Fortran', 'Programming', 'Language']) , deliml = SK_"""" )
    call disp%skip()

end program example