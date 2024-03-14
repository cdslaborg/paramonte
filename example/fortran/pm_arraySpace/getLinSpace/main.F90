program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_kind, only: RK32, RK64, RK128, CK32, CK64, CK128
    use pm_arraySpace, only: getLinSpace 

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Generate real linspace.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!32-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace(0._RK32, 10._RK32, 4_IK)")
    call disp%show( getLinSpace(0._RK32, 10._RK32, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(0._RK32, 10._RK32, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(0._RK32, 10._RK32, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RK32, 10._RK32, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RK32, 10._RK32, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RK32, 10._RK32, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RK32, 10._RK32, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace(10._RK32, 0._RK32, 4_IK)")
    call disp%show( getLinSpace(10._RK32, 0._RK32, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(10._RK32, 0._RK32, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(10._RK32, 0._RK32, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RK32, 0._RK32, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RK32, 0._RK32, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RK32, 0._RK32, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RK32, 0._RK32, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!64-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace(0._RK64, 10._RK64, 4_IK)")
    call disp%show( getLinSpace(0._RK64, 10._RK64, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(0._RK64, 10._RK64, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(0._RK64, 10._RK64, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RK64, 10._RK64, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RK64, 10._RK64, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RK64, 10._RK64, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RK64, 10._RK64, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace(10._RK64, 0._RK64, 4_IK)")
    call disp%show( getLinSpace(10._RK64, 0._RK64, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(10._RK64, 0._RK64, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(10._RK64, 0._RK64, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RK64, 0._RK64, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RK64, 0._RK64, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RK64, 0._RK64, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RK64, 0._RK64, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!128-bit real")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace(0._RK128, 10._RK128, 4_IK)")
    call disp%show( getLinSpace(0._RK128, 10._RK128, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(0._RK128, 10._RK128, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(0._RK128, 10._RK128, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RK128, 10._RK128, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RK128, 10._RK128, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RK128, 10._RK128, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RK128, 10._RK128, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace(10._RK128, 0._RK128, 4_IK)")
    call disp%show( getLinSpace(10._RK128, 0._RK128, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(10._RK128, 0._RK128, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(10._RK128, 0._RK128, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RK128, 0._RK128, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RK128, 0._RK128, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RK128, 0._RK128, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RK128, 0._RK128, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip


    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Generate complex linspace.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip


    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!32-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace((0._CK32,0._CK32), (10._CK32,10._CK32), 4_IK)")
    call disp%show( getLinSpace((0._CK32,0._CK32), (10._CK32,10._CK32), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((0._CK32,0._CK32), (10._CK32,10._CK32), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((0._CK32,0._CK32), (10._CK32,10._CK32), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CK32,0._CK32), (10._CK32,10._CK32), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CK32,0._CK32), (10._CK32,10._CK32), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CK32,0._CK32), (10._CK32,10._CK32), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CK32,0._CK32), (10._CK32,10._CK32), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace((10._CK32,-10._CK32), (0._CK32,0._CK32), 4_IK)")
    call disp%show( getLinSpace((10._CK32,-10._CK32), (0._CK32,0._CK32), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((10._CK32,-10._CK32), (0._CK32,0._CK32), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((10._CK32,-10._CK32), (0._CK32,0._CK32), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CK32,-10._CK32), (0._CK32,0._CK32), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CK32,-10._CK32), (0._CK32,0._CK32), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CK32,-10._CK32), (0._CK32,0._CK32), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CK32,-10._CK32), (0._CK32,0._CK32), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace((0._CK64,0._CK64), (10._CK64,10._CK64), 4_IK)")
    call disp%show( getLinSpace((0._CK64,0._CK64), (10._CK64,10._CK64), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((0._CK64,0._CK64), (10._CK64,10._CK64), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((0._CK64,0._CK64), (10._CK64,10._CK64), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CK64,0._CK64), (10._CK64,10._CK64), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CK64,0._CK64), (10._CK64,10._CK64), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CK64,0._CK64), (10._CK64,10._CK64), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CK64,0._CK64), (10._CK64,10._CK64), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace((10._CK64,-10._CK64), (0._CK64,0._CK64), 4_IK)")
    call disp%show( getLinSpace((10._CK64,-10._CK64), (0._CK64,0._CK64), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((10._CK64,-10._CK64), (0._CK64,0._CK64), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((10._CK64,-10._CK64), (0._CK64,0._CK64), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CK64,-10._CK64), (0._CK64,0._CK64), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CK64,-10._CK64), (0._CK64,0._CK64), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CK64,-10._CK64), (0._CK64,0._CK64), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CK64,-10._CK64), (0._CK64,0._CK64), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace((0._CK128,0._CK128), (10._CK128,10._CK128), 4_IK)")
    call disp%show( getLinSpace((0._CK128,0._CK128), (10._CK128,10._CK128), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((0._CK128,0._CK128), (10._CK128,10._CK128), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((0._CK128,0._CK128), (10._CK128,10._CK128), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CK128,0._CK128), (10._CK128,10._CK128), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CK128,0._CK128), (10._CK128,10._CK128), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CK128,0._CK128), (10._CK128,10._CK128), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CK128,0._CK128), (10._CK128,10._CK128), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace((10._CK128,-10._CK128), (0._CK128,0._CK128), 4_IK)")
    call disp%show( getLinSpace((10._CK128,-10._CK128), (0._CK128,0._CK128), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((10._CK128,-10._CK128), (0._CK128,0._CK128), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((10._CK128,-10._CK128), (0._CK128,0._CK128), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CK128,-10._CK128), (0._CK128,0._CK128), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CK128,-10._CK128), (0._CK128,0._CK128), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CK128,-10._CK128), (0._CK128,0._CK128), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CK128,-10._CK128), (0._CK128,0._CK128), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

end program example