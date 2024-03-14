program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: getFileUnit

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Get the unit of a file that is already opened.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("open(unit = 11, file = 'temp.tmp')")
                    open(unit = 11, file = 'temp.tmp')
    call disp%show("getFileUnit('temp.tmp')")
    call disp%show( getFileUnit('temp.tmp') )
    call disp%show("close(unit = 11)")
                    close(unit = 11)
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Get the units of multiple files.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("open(unit = 11, file = 'temp.tmp')")
                    open(unit = 11, file = 'temp.tmp')
    call disp%show("getFileUnit(['temp.tmp', 'tmp.temp']) ! The second file is not connected., yielding a unit of `-1`.")
    call disp%show( getFileUnit(['temp.tmp', 'tmp.temp']) )
    call disp%show("close(unit = 11)")
                    close(unit = 11)
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Get a new unique unit number for use in an `open()` statement.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getFileUnit()")
    call disp%show( getFileUnit() )
    call disp%skip

end program example