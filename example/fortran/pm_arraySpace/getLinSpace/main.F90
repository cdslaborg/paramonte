program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_kind, only: RKS, RKD, RKH, CKS, CKD, CKH
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
    call disp%show("getLinSpace(0._RKS, 10._RKS, 4_IK)")
    call disp%show( getLinSpace(0._RKS, 10._RKS, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(0._RKS, 10._RKS, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(0._RKS, 10._RKS, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RKS, 10._RKS, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RKS, 10._RKS, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RKS, 10._RKS, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RKS, 10._RKS, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace(10._RKS, 0._RKS, 4_IK)")
    call disp%show( getLinSpace(10._RKS, 0._RKS, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(10._RKS, 0._RKS, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(10._RKS, 0._RKS, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RKS, 0._RKS, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RKS, 0._RKS, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RKS, 0._RKS, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RKS, 0._RKS, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!64-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace(0._RKD, 10._RKD, 4_IK)")
    call disp%show( getLinSpace(0._RKD, 10._RKD, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(0._RKD, 10._RKD, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(0._RKD, 10._RKD, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RKD, 10._RKD, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RKD, 10._RKD, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RKD, 10._RKD, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RKD, 10._RKD, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace(10._RKD, 0._RKD, 4_IK)")
    call disp%show( getLinSpace(10._RKD, 0._RKD, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(10._RKD, 0._RKD, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(10._RKD, 0._RKD, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RKD, 0._RKD, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RKD, 0._RKD, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RKD, 0._RKD, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RKD, 0._RKD, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!128-bit real")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace(0._RKH, 10._RKH, 4_IK)")
    call disp%show( getLinSpace(0._RKH, 10._RKH, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(0._RKH, 10._RKH, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(0._RKH, 10._RKH, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RKH, 10._RKH, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RKH, 10._RKH, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(0._RKH, 10._RKH, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(0._RKH, 10._RKH, 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace(10._RKH, 0._RKH, 4_IK)")
    call disp%show( getLinSpace(10._RKH, 0._RKH, 4_IK) )
    call disp%skip
    call disp%show("getLinSpace(10._RKH, 0._RKH, 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace(10._RKH, 0._RKH, 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RKH, 0._RKH, 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RKH, 0._RKH, 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace(10._RKH, 0._RKH, 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace(10._RKH, 0._RKH, 4_IK, fopen = .true._LK, lopen = .true._LK) )
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
    call disp%show("getLinSpace((0._CKS,0._CKS), (10._CKS,10._CKS), 4_IK)")
    call disp%show( getLinSpace((0._CKS,0._CKS), (10._CKS,10._CKS), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((0._CKS,0._CKS), (10._CKS,10._CKS), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((0._CKS,0._CKS), (10._CKS,10._CKS), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CKS,0._CKS), (10._CKS,10._CKS), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CKS,0._CKS), (10._CKS,10._CKS), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CKS,0._CKS), (10._CKS,10._CKS), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CKS,0._CKS), (10._CKS,10._CKS), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace((10._CKS,-10._CKS), (0._CKS,0._CKS), 4_IK)")
    call disp%show( getLinSpace((10._CKS,-10._CKS), (0._CKS,0._CKS), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((10._CKS,-10._CKS), (0._CKS,0._CKS), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((10._CKS,-10._CKS), (0._CKS,0._CKS), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CKS,-10._CKS), (0._CKS,0._CKS), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CKS,-10._CKS), (0._CKS,0._CKS), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CKS,-10._CKS), (0._CKS,0._CKS), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CKS,-10._CKS), (0._CKS,0._CKS), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace((0._CKD,0._CKD), (10._CKD,10._CKD), 4_IK)")
    call disp%show( getLinSpace((0._CKD,0._CKD), (10._CKD,10._CKD), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((0._CKD,0._CKD), (10._CKD,10._CKD), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((0._CKD,0._CKD), (10._CKD,10._CKD), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CKD,0._CKD), (10._CKD,10._CKD), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CKD,0._CKD), (10._CKD,10._CKD), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CKD,0._CKD), (10._CKD,10._CKD), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CKD,0._CKD), (10._CKD,10._CKD), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace((10._CKD,-10._CKD), (0._CKD,0._CKD), 4_IK)")
    call disp%show( getLinSpace((10._CKD,-10._CKD), (0._CKD,0._CKD), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((10._CKD,-10._CKD), (0._CKD,0._CKD), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((10._CKD,-10._CKD), (0._CKD,0._CKD), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CKD,-10._CKD), (0._CKD,0._CKD), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CKD,-10._CKD), (0._CKD,0._CKD), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CKD,-10._CKD), (0._CKD,0._CKD), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CKD,-10._CKD), (0._CKD,0._CKD), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLinSpace((0._CKH,0._CKH), (10._CKH,10._CKH), 4_IK)")
    call disp%show( getLinSpace((0._CKH,0._CKH), (10._CKH,10._CKH), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((0._CKH,0._CKH), (10._CKH,10._CKH), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((0._CKH,0._CKH), (10._CKH,10._CKH), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CKH,0._CKH), (10._CKH,10._CKH), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CKH,0._CKH), (10._CKH,10._CKH), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((0._CKH,0._CKH), (10._CKH,10._CKH), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((0._CKH,0._CKH), (10._CKH,10._CKH), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLinSpace((10._CKH,-10._CKH), (0._CKH,0._CKH), 4_IK)")
    call disp%show( getLinSpace((10._CKH,-10._CKH), (0._CKH,0._CKH), 4_IK) )
    call disp%skip
    call disp%show("getLinSpace((10._CKH,-10._CKH), (0._CKH,0._CKH), 4_IK, fopen = .true._LK)")
    call disp%show( getLinSpace((10._CKH,-10._CKH), (0._CKH,0._CKH), 4_IK, fopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CKH,-10._CKH), (0._CKH,0._CKH), 4_IK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CKH,-10._CKH), (0._CKH,0._CKH), 4_IK, lopen = .true._LK) )
    call disp%skip
    call disp%show("getLinSpace((10._CKH,-10._CKH), (0._CKH,0._CKH), 4_IK, fopen = .true._LK, lopen = .true._LK)")
    call disp%show( getLinSpace((10._CKH,-10._CKH), (0._CKH,0._CKH), 4_IK, fopen = .true._LK, lopen = .true._LK) )
    call disp%skip

end program example