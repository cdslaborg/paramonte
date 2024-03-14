program example

    use pm_kind, only: LK, SK
    use pm_io, only: display_type
    use pm_kind, only: IK, RK32, RK64, RK128, CK32, CK64, CK128
    use pm_arraySpace, only: getLogSpace 

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")


    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("Generate real logspace.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")

    call disp%skip
    call disp%show("!%%%%%%%%%%")
    call disp%show("32-bit real")
    call disp%show("!%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace(0._RK32, 4._RK32, 4_IK)")
    call disp%show( getLogSpace(0._RK32, 4._RK32, 4_IK) )
    call disp%show("getLogSpace(0._RK32, 4._RK32, 4_IK, base = 10._RK32)")
    call disp%show( getLogSpace(0._RK32, 4._RK32, 4_IK, base = 10._RK32) )
    call disp%show("getLogSpace(0._RK32, 4._RK32, 4_IK, fopen = .true._LK, base = 10._RK32)")
    call disp%show( getLogSpace(0._RK32, 4._RK32, 4_IK, fopen = .true._LK, base = 10._RK32) )
    call disp%show("getLogSpace(0._RK32, 4._RK32, 4_IK, lopen = .true._LK, base = 10._RK32)")
    call disp%show( getLogSpace(0._RK32, 4._RK32, 4_IK, lopen = .true._LK, base = 10._RK32) )
    call disp%show("getLogSpace(0._RK32, 4._RK32, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK32)")
    call disp%show( getLogSpace(0._RK32, 4._RK32, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK32) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace(4._RK32, 0._RK32, 4_IK)")
    call disp%show( getLogSpace(4._RK32, 0._RK32, 4_IK) )
    call disp%show("getLogSpace(4._RK32, 0._RK32, 4_IK, base = 10._RK32)")
    call disp%show( getLogSpace(4._RK32, 0._RK32, 4_IK, base = 10._RK32) )
    call disp%show("getLogSpace(4._RK32, 0._RK32, 4_IK, fopen = .true._LK, base = 10._RK32)")
    call disp%show( getLogSpace(4._RK32, 0._RK32, 4_IK, fopen = .true._LK, base = 10._RK32) )
    call disp%show("getLogSpace(4._RK32, 0._RK32, 4_IK, lopen = .true._LK, base = 10._RK32)")
    call disp%show( getLogSpace(4._RK32, 0._RK32, 4_IK, lopen = .true._LK, base = 10._RK32) )
    call disp%show("getLogSpace(4._RK32, 0._RK32, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK32)")
    call disp%show( getLogSpace(4._RK32, 0._RK32, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK32) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%")
    call disp%show("64-bit real")
    call disp%show("!%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace(0._RK64, 4._RK64, 4_IK)")
    call disp%show( getLogSpace(0._RK64, 4._RK64, 4_IK) )
    call disp%show("getLogSpace(0._RK64, 4._RK64, 4_IK, base = 10._RK64)")
    call disp%show( getLogSpace(0._RK64, 4._RK64, 4_IK, base = 10._RK64) )
    call disp%show("getLogSpace(0._RK64, 4._RK64, 4_IK, fopen = .true._LK, base = 10._RK64)")
    call disp%show( getLogSpace(0._RK64, 4._RK64, 4_IK, fopen = .true._LK, base = 10._RK64) )
    call disp%show("getLogSpace(0._RK64, 4._RK64, 4_IK, lopen = .true._LK, base = 10._RK64)")
    call disp%show( getLogSpace(0._RK64, 4._RK64, 4_IK, lopen = .true._LK, base = 10._RK64) )
    call disp%show("getLogSpace(0._RK64, 4._RK64, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK64)")
    call disp%show( getLogSpace(0._RK64, 4._RK64, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK64) )
    call disp%skip
    
    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace(4._RK64, 0._RK64, 4_IK)")
    call disp%show( getLogSpace(4._RK64, 0._RK64, 4_IK) )
    call disp%show("getLogSpace(4._RK64, 0._RK64, 4_IK, base = 10._RK64)")
    call disp%show( getLogSpace(4._RK64, 0._RK64, 4_IK, base = 10._RK64) )
    call disp%show("getLogSpace(4._RK64, 0._RK64, 4_IK, fopen = .true._LK, base = 10._RK64)")
    call disp%show( getLogSpace(4._RK64, 0._RK64, 4_IK, fopen = .true._LK, base = 10._RK64) )
    call disp%show("getLogSpace(4._RK64, 0._RK64, 4_IK, lopen = .true._LK, base = 10._RK64)")
    call disp%show( getLogSpace(4._RK64, 0._RK64, 4_IK, lopen = .true._LK, base = 10._RK64) )
    call disp%show("getLogSpace(4._RK64, 0._RK64, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK64)")
    call disp%show( getLogSpace(4._RK64, 0._RK64, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK64) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("128-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace(0._RK128, 4._RK128, 4_IK)")
    call disp%show( getLogSpace(0._RK128, 4._RK128, 4_IK) )
    call disp%show("getLogSpace(0._RK128, 4._RK128, 4_IK, base = 10._RK128)")
    call disp%show( getLogSpace(0._RK128, 4._RK128, 4_IK, base = 10._RK128) )
    call disp%show("getLogSpace(0._RK128, 4._RK128, 4_IK, fopen = .true._LK, base = 10._RK128)")
    call disp%show( getLogSpace(0._RK128, 4._RK128, 4_IK, fopen = .true._LK, base = 10._RK128) )
    call disp%show("getLogSpace(0._RK128, 4._RK128, 4_IK, lopen = .true._LK, base = 10._RK128)")
    call disp%show( getLogSpace(0._RK128, 4._RK128, 4_IK, lopen = .true._LK, base = 10._RK128) )
    call disp%show("getLogSpace(0._RK128, 4._RK128, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK128)")
    call disp%show( getLogSpace(0._RK128, 4._RK128, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK128) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace(4._RK128, 0._RK128, 4_IK)")
    call disp%show( getLogSpace(4._RK128, 0._RK128, 4_IK) )
    call disp%show("getLogSpace(4._RK128, 0._RK128, 4_IK, base = 10._RK128)")
    call disp%show( getLogSpace(4._RK128, 0._RK128, 4_IK, base = 10._RK128) )
    call disp%show("getLogSpace(4._RK128, 0._RK128, 4_IK, fopen = .true._LK, base = 10._RK128)")
    call disp%show( getLogSpace(4._RK128, 0._RK128, 4_IK, fopen = .true._LK, base = 10._RK128) )
    call disp%show("getLogSpace(4._RK128, 0._RK128, 4_IK, lopen = .true._LK, base = 10._RK128)")
    call disp%show( getLogSpace(4._RK128, 0._RK128, 4_IK, lopen = .true._LK, base = 10._RK128) )
    call disp%show("getLogSpace(4._RK128, 0._RK128, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK128)")
    call disp%show( getLogSpace(4._RK128, 0._RK128, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RK128) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("Generate complex logspace.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("32-bit complex")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace((1._CK32,1._CK32), (4._CK32,4._CK32), 4_IK)")
    call disp%show( getLogSpace((1._CK32,1._CK32), (4._CK32,4._CK32), 4_IK) )
    call disp%show("getLogSpace((1._CK32,1._CK32), (4._CK32,4._CK32), 4_IK, base = 10._CK32)")
    call disp%show( getLogSpace((1._CK32,1._CK32), (4._CK32,4._CK32), 4_IK, base = 10._CK32) )
    call disp%show("getLogSpace((0._CK32,0._CK32), (4._CK32,4._CK32), 4_IK, fopen = .true._LK, base = 10._CK32)")
    call disp%show( getLogSpace((0._CK32,0._CK32), (4._CK32,4._CK32), 4_IK, fopen = .true._LK, base = 10._CK32) )
    call disp%show("getLogSpace((0._CK32,0._CK32), (4._CK32,4._CK32), 4_IK, lopen = .true._LK, base = 10._CK32)")
    call disp%show( getLogSpace((0._CK32,0._CK32), (4._CK32,4._CK32), 4_IK, lopen = .true._LK, base = 10._CK32) )
    call disp%show("getLogSpace((0._CK32,0._CK32), (4._CK32,4._CK32), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK32)")
    call disp%show( getLogSpace((0._CK32,0._CK32), (4._CK32,4._CK32), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK32) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK)")
    call disp%show( getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK) )
    call disp%show("getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK, base = 10._CK32)")
    call disp%show( getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK, base = 10._CK32) )
    call disp%show("getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK, fopen = .true._LK, base = 10._CK32)")
    call disp%show( getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK, fopen = .true._LK, base = 10._CK32) )
    call disp%show("getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK, lopen = .true._LK, base = 10._CK32)")
    call disp%show( getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK, lopen = .true._LK, base = 10._CK32) )
    call disp%show("getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK32)")
    call disp%show( getLogSpace((4._CK32,-4._CK32), (-4._CK32,+4._CK32), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK32) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("64-bit complex")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK)")
    call disp%show( getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK) )
    call disp%show("getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK, base = 10._CK64)")
    call disp%show( getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK, base = 10._CK64) )
    call disp%show("getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK, fopen = .true._LK, base = 10._CK64)")
    call disp%show( getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK, fopen = .true._LK, base = 10._CK64) )
    call disp%show("getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK, lopen = .true._LK, base = 10._CK64)")
    call disp%show( getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK, lopen = .true._LK, base = 10._CK64) )
    call disp%show("getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK64)")
    call disp%show( getLogSpace((0._CK64,0._CK64), (4._CK64,4._CK64), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK64) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK)")
    call disp%show( getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK) )
    call disp%show("getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK, base = 10._CK64)")
    call disp%show( getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK, base = 10._CK64) )
    call disp%show("getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK, fopen = .true._LK, base = 10._CK64)")
    call disp%show( getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK, fopen = .true._LK, base = 10._CK64) )
    call disp%show("getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK, lopen = .true._LK, base = 10._CK64)")
    call disp%show( getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK, lopen = .true._LK, base = 10._CK64) )
    call disp%show("getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK64)")
    call disp%show( getLogSpace((4._CK64,-4._CK64), (-4._CK64,+4._CK64), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK64) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK)")
    call disp%show( getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK) )
    call disp%show("getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK, base = 10._CK128)")
    call disp%show( getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK, base = 10._CK128) )
    call disp%show("getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK, fopen = .true._LK, base = 10._CK128)")
    call disp%show( getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK, fopen = .true._LK, base = 10._CK128) )
    call disp%show("getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK, lopen = .true._LK, base = 10._CK128)")
    call disp%show( getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK, lopen = .true._LK, base = 10._CK128) )
    call disp%show("getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK128)")
    call disp%show( getLogSpace((0._CK128,0._CK128), (4._CK128,4._CK128), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK128) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK)")
    call disp%show( getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK) )
    call disp%show("getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK, base = 10._CK128)")
    call disp%show( getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK, base = 10._CK128) )
    call disp%show("getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK, fopen = .true._LK, base = 10._CK128)")
    call disp%show( getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK, fopen = .true._LK, base = 10._CK128) )
    call disp%show("getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK, lopen = .true._LK, base = 10._CK128)")
    call disp%show( getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK, lopen = .true._LK, base = 10._CK128) )
    call disp%show("getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK128)")
    call disp%show( getLogSpace((4._CK128,-4._CK128), (-4._CK128,+4._CK128), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CK128) )
    call disp%skip

end program example