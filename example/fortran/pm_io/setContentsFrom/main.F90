program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: setContentsFrom
    use pm_arrayCenter, only: getCentered

    implicit none

    integer(IK)                     :: iostat
    character(255, SK)              :: iomsg
    character(  :, SK), allocatable :: contents

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("call setContentsFrom('main.F90', contents, iostat = iostat, iomsg = iomsg)")
                    call setContentsFrom('main.F90', contents, iostat = iostat, iomsg = iomsg)
    call disp%show( getCentered(SK_" contents ", size = 132_IK, fill = SK_"_") )
    call disp%show( contents )
    call disp%show( getCentered(SK_" end of contents ", size = 132_IK, fill = SK_"_") )
    call disp%show("iostat")
    call disp%show( iostat )
    call disp%show("trim(iomsg)")
    call disp%show( trim(iomsg) , deliml = SK_"""" )
    call disp%skip

end program example