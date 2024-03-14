program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: setFileClosed

    implicit none

    integer(IK)         :: iostat
    character(132,SK)   :: iomsg

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Control runtime IO error with optional error status code and message.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    iomsg = repeat(" ", len(iomsg))

    call disp%skip
    call disp%show("call setFileClosed(unit = 43, iostat = iostat, iomsg = iomsg)")
                    call setFileClosed(unit = 43, iostat = iostat, iomsg = iomsg)
    call disp%show("iostat")
    call disp%show( iostat )
    call disp%show("trim(iomsg)")
    call disp%show( trim(iomsg) , deliml = """")
    call disp%skip

    iomsg = repeat(" ", len(iomsg))

    call disp%skip
    call disp%show("call setFileClosed(unit = 43, del = .true._LK, iostat = iostat, iomsg = iomsg)")
                    call setFileClosed(unit = 43, del = .true._LK, iostat = iostat, iomsg = iomsg)
    call disp%show("iostat")
    call disp%show( iostat )
    call disp%show("trim(iomsg)")
    call disp%show( trim(iomsg) , deliml = """")
    call disp%skip

    iomsg = repeat(" ", len(iomsg))

    call disp%skip
    call disp%show("open(unit = 34, status = 'scratch')")
                    open(unit = 34, status = 'scratch')
    call disp%show("call setFileClosed(unit = 34, del = .true._LK, iostat = iostat, iomsg = iomsg)")
                    call setFileClosed(unit = 34, del = .true._LK, iostat = iostat, iomsg = iomsg)
    call disp%show("iostat")
    call disp%show( iostat )
    call disp%show("trim(iomsg)")
    call disp%show( trim(iomsg) , deliml = """")
    call disp%skip

    iomsg = repeat(" ", len(iomsg))

    call disp%skip
    call disp%show("call setFileClosed(unit = -34, iostat = iostat)")
                    call setFileClosed(unit = -34, iostat = iostat)
    call disp%show("iostat")
    call disp%show( iostat )
    call disp%show("trim(iomsg)")
    call disp%show( trim(iomsg) , deliml = SK_"""" )
    call disp%skip

end program example