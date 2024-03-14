program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: setRecordFrom

    implicit none

    integer(IK)                     :: ub
    integer(IK)                     :: unit
    integer(IK)                     :: iostat
    character(255, SK)              :: iomsg
    character(  :, SK), allocatable :: record

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Read the record-oriented (sequential access) `main.F90` file, line by line to the end.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    iomsg = repeat(" ", len(iomsg))

    open(newunit = unit, file = "main.F90", status = "old")

    do
        call disp%skip
        call disp%show("call setRecordFrom(unit, record, iostat = iostat, iomsg = iomsg)")
                        call setRecordFrom(unit, record, iostat = iostat, iomsg = iomsg)
        if (iostat == 0_IK) then
            call disp%show("record")
            call disp%show( record , deliml = SK_"""" )
        else
            call disp%show("iostat")
            call disp%show( iostat )
            call disp%show("trim(adjustl(iomsg))")
            call disp%show( trim(adjustl(iomsg)) , deliml = SK_"""" )
            exit
        end if
        call disp%skip
    end do

    close(unit)

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Read a record-oriented (sequential access) `main.F90` file, line by line to the end, faster without redundant reallocations and copies.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    iomsg = repeat(" ", len(iomsg))

    call disp%show('open(newunit = unit, file = "main.F90", status = "old")')
                    open(newunit = unit, file = "main.F90", status = "old")
    do
        call disp%skip
        call disp%show("call setRecordFrom(unit, record, iostat = iostat, iomsg = iomsg, ub = ub)")
                        call setRecordFrom(unit, record, iostat = iostat, iomsg = iomsg, ub = ub)
        if (iostat == 0_IK) then
            call disp%show("record(1:ub)")
            call disp%show( record(1:ub) , deliml = SK_"""" )
        else
            call disp%show("iostat")
            call disp%show( iostat )
            call disp%show("trim(adjustl(iomsg))")
            call disp%show( trim(adjustl(iomsg)) , deliml = SK_"""" )
            exit
        end if
        call disp%skip
    end do

    close(unit)

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Read a record-oriented (sequential access) `main.F90` file line by line to the end in a single string, each line separated by linefeed (see the end of contents for the source code).")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block
        use iso_fortran_env, only: iostat_end
        use pm_arrayCenter, only: getCentered
        integer(IK) :: lb, ub
        open(newunit = unit, file = "main.F90", status = "old")
        iomsg = repeat(" ", len(iomsg))
        lb = 1_IK
        do
            call setRecordFrom(unit, record, iostat = iostat, iomsg = iomsg, lb = lb, ub = ub, linefed = .true._LK)
            if (iostat == 0_IK) then
                lb = ub + 1_IK
                cycle
            elseif (iostat == iostat_end) then
                exit
            else
                error stop SK_"Unknown IO error occurred: "//trim(iomsg)
            end if
        end do
        call disp%show( getCentered(SK_" records ", size = 132_IK, fill = SK_"_") )
        call disp%show( record(1:ub) )
        call disp%show( getCentered(SK_" end of records ", size = 132_IK, fill = SK_"_") )
        close(unit)
    end block

end program example