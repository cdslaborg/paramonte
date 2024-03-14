program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: getRecordFrom
    use pm_io, only: display_type
    use iso_fortran_env, only: iostat_end

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

    call disp%show('open(newunit = unit, file = "main.F90", status = "old")')
                    open(newunit = unit, file = "main.F90", status = "old")
    do
        call disp%skip
        call disp%show("record = getRecordFrom(unit, iostat = iostat)")
                        record = getRecordFrom(unit, iostat = iostat)
        if (iostat == 0_IK) then
            call disp%show("record")
            call disp%show( record , deliml = SK_"""" )
        else
            call disp%show("[iostat, iostat_end]")
            call disp%show( [iostat, iostat_end] )
            exit
        end if
        call disp%skip
    end do

    close(unit)

end program example