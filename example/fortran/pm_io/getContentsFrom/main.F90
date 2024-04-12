program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: getContentsFrom
    use pm_arrayCenter, only: getCentered

    implicit none

    integer :: unit
    character(:, SK), allocatable :: contents

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("contents = getContentsFrom('main.F90')")
                    contents = getContentsFrom('main.F90')
    call disp%show( getCentered(SK_" contents ", size = 132_IK, fill = SK_"_") )
    call disp%show( contents )
    call disp%show( getCentered(SK_" end of contents ", size = 132_IK, fill = SK_"_") )
    call disp%skip

    ! Intel ifort Windows bug?
    ! @file(pm_io@routines.inc.F90)@line(235)@pm_io@setContentsFrom(): A non-advancing READ immediately following a non-advancing WRITE on the same unit number is not allowed, unit -129.
#if !__INTEL_COMPILER || !_WIN32
    call disp%skip
    call disp%show("open(newunit = unit, file = 'main.F90')")
                    open(newunit = unit, file = 'main.F90')
    call disp%show("contents = getContentsFrom(unit)")
                    contents = getContentsFrom(unit)
    call disp%show( getCentered(SK_" contents ", size = 132_IK, fill = SK_"_") )
    call disp%show( contents )
    call disp%show( getCentered(SK_" end of contents ", size = 132_IK, fill = SK_"_") )
    call disp%skip
#endif

end program example