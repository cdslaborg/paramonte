program example

    use pm_str, only: NLC
    use pm_val2str, only: getStr
    use pm_kind, only: SK, IK, LK
    use pm_arrayCenter, only: getCentered
    use pm_io, only: display_type
    use pm_io, only: getContentsFrom
    use pm_io, only: setContentsTo

    implicit none

    integer(IK)                     :: iostat
    character(255, SK)              :: iomsg
    character(  :, SK), allocatable :: contents

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("contents = repeat(SK_'""a"", '//getStr([1, 2, 3])//NLC, 2_IK)")
                    contents = repeat(SK_'"a", '//getStr([1, 2, 3])//NLC, 2_IK)
    call disp%show("call setContentsTo('temp.temp', contents)")
                    call setContentsTo('temp.temp', contents)
    call disp%show("call setContentsTo('temp.temp', contents, iostat = iostat, iomsg = iomsg) ! optionally catch io errors.")
                    call setContentsTo('temp.temp', contents, iostat = iostat, iomsg = iomsg)
    call disp%show("iostat")
    call disp%show( iostat )
    call disp%show("trim(iomsg)")
    call disp%show( trim(iomsg) , deliml = SK_"""" )
    call disp%show("getContentsFrom('temp.temp')")
    call disp%show( getContentsFrom('temp.temp') )
    call disp%skip

end program example