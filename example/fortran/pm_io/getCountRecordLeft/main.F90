program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: getCountRecordLeft

    implicit none

    integer(IK) :: unit
    character(255, SK) :: record
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Count the number of records in the record-oriented (sequential access) `main.F90` file.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("open(newunit = unit, file = 'main.F90', status = 'old', action = 'read')")
                    open(newunit = unit, file = 'main.F90', status = 'old', action = 'read')
    call disp%show("getCountRecordLeft(unit)")
    call disp%show( getCountRecordLeft(unit) )
    call disp%show("close(unit)")
                    close(unit)
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Count the number of records in the `main.F90` file starting from the current position, then reset the position to where it was.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("open(newunit = unit, file = 'main.F90', status = 'old', action = 'read')")
                    open(newunit = unit, file = 'main.F90', status = 'old', action = 'read')
    call disp%show("read(unit,*); read(unit,*) ! move the position to the beginning of the third record.")
                    read(unit,*); read(unit,*)
    call disp%show("read(unit,'(a)') record; backspace(unit) ! read the current line as a reference to compare.")
                    read(unit,'(a)') record; backspace(unit)
    call disp%show("trim(record)")
    call disp%show( trim(record) , deliml = SK_"""" )
    call disp%show("getCountRecordLeft(unit, reset = .true._LK)")
    call disp%show( getCountRecordLeft(unit, reset = .true._LK) )
    call disp%show("read(unit,'(a)') record ! read the current line again to compare with the value before calling getCountRecordLeft().")
                    read(unit,'(a)') record
    call disp%show("trim(record)")
    call disp%show( trim(record) , deliml = SK_"""" )
    call disp%show("close(unit)")
                    close(unit)
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Exclude certain records that match the user-specified behavior via isCountable.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("open(newunit = unit, file = 'main.F90', status = 'old', action = 'read')")
                    open(newunit = unit, file = 'main.F90', status = 'old', action = 'read')
    call disp%show("getCountRecordLeft(unit, isCountable = isCountable)")
    call disp%show( getCountRecordLeft(unit, isCountable = isCountable) )
    call disp%show("close(unit)")
                    close(unit)
    call disp%skip

contains

    function isCountable(record) result(countable)
        character(*, SK), intent(in)    :: record
        logical(LK)                     :: countable
        countable = index(record, SK_"call disp%skip", kind = IK) == 0_IK
    end function

end program example