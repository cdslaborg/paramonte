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

    call disp%skip
    call disp%show("open(newunit = unit, file = 'main.F90')")
                    open(newunit = unit, file = 'main.F90')
    call disp%show("contents = getContentsFrom(unit)")
                    contents = getContentsFrom(unit)
    call disp%show( getCentered(SK_" contents ", size = 132_IK, fill = SK_"_") )
    call disp%show( contents )
    call disp%show( getCentered(SK_" end of contents ", size = 132_IK, fill = SK_"_") )
    call disp%skip

end program example