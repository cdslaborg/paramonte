program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: getCountRecord

    implicit none

    integer(IK)                     :: iostat
    character(255, SK)              :: iomsg

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Count the number of records in the record-oriented (sequential access) `main.F90` file.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    iomsg = repeat(" ", len(iomsg))

    call disp%skip
    call disp%show("getCountRecord(file = SK_'main.F90')")
    call disp%show( getCountRecord(file = SK_'main.F90') )
    call disp%skip

    call disp%skip
    call disp%show("getCountRecord(file = SK_'main.F90', iostat = iostat, iomsg = iomsg)")
    call disp%show( getCountRecord(file = SK_'main.F90', iostat = iostat, iomsg = iomsg) )
    call disp%show("iostat")
    call disp%show( iostat )
    call disp%show("trim(adjustl(iomsg))")
    call disp%show( trim(adjustl(iomsg)) , deliml = SK_"""" )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Exclude certain records that match the user-specified behavior via isCountable.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getCountRecord(file = SK_'main.F90', isCountable = isCountable)")
    call disp%show( getCountRecord(file = SK_'main.F90', isCountable = isCountable) )
    call disp%skip

contains

    function isCountable(record) result(countable)
        character(*, SK), intent(in)    :: record
        logical(LK)                     :: countable
        countable = index(record, SK_"call disp%skip", kind = IK) == 0_IK
    end function

end program example