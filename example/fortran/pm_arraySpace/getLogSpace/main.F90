program example

    use pm_kind, only: LK, SK
    use pm_io, only: display_type
    use pm_kind, only: IK, RKS, RKD, RKH, CKS, CKD, CKH
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
    call disp%show("getLogSpace(0._RKS, 4._RKS, 4_IK)")
    call disp%show( getLogSpace(0._RKS, 4._RKS, 4_IK) )
    call disp%show("getLogSpace(0._RKS, 4._RKS, 4_IK, base = 10._RKS)")
    call disp%show( getLogSpace(0._RKS, 4._RKS, 4_IK, base = 10._RKS) )
    call disp%show("getLogSpace(0._RKS, 4._RKS, 4_IK, fopen = .true._LK, base = 10._RKS)")
    call disp%show( getLogSpace(0._RKS, 4._RKS, 4_IK, fopen = .true._LK, base = 10._RKS) )
    call disp%show("getLogSpace(0._RKS, 4._RKS, 4_IK, lopen = .true._LK, base = 10._RKS)")
    call disp%show( getLogSpace(0._RKS, 4._RKS, 4_IK, lopen = .true._LK, base = 10._RKS) )
    call disp%show("getLogSpace(0._RKS, 4._RKS, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKS)")
    call disp%show( getLogSpace(0._RKS, 4._RKS, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKS) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace(4._RKS, 0._RKS, 4_IK)")
    call disp%show( getLogSpace(4._RKS, 0._RKS, 4_IK) )
    call disp%show("getLogSpace(4._RKS, 0._RKS, 4_IK, base = 10._RKS)")
    call disp%show( getLogSpace(4._RKS, 0._RKS, 4_IK, base = 10._RKS) )
    call disp%show("getLogSpace(4._RKS, 0._RKS, 4_IK, fopen = .true._LK, base = 10._RKS)")
    call disp%show( getLogSpace(4._RKS, 0._RKS, 4_IK, fopen = .true._LK, base = 10._RKS) )
    call disp%show("getLogSpace(4._RKS, 0._RKS, 4_IK, lopen = .true._LK, base = 10._RKS)")
    call disp%show( getLogSpace(4._RKS, 0._RKS, 4_IK, lopen = .true._LK, base = 10._RKS) )
    call disp%show("getLogSpace(4._RKS, 0._RKS, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKS)")
    call disp%show( getLogSpace(4._RKS, 0._RKS, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKS) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%")
    call disp%show("64-bit real")
    call disp%show("!%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace(0._RKD, 4._RKD, 4_IK)")
    call disp%show( getLogSpace(0._RKD, 4._RKD, 4_IK) )
    call disp%show("getLogSpace(0._RKD, 4._RKD, 4_IK, base = 10._RKD)")
    call disp%show( getLogSpace(0._RKD, 4._RKD, 4_IK, base = 10._RKD) )
    call disp%show("getLogSpace(0._RKD, 4._RKD, 4_IK, fopen = .true._LK, base = 10._RKD)")
    call disp%show( getLogSpace(0._RKD, 4._RKD, 4_IK, fopen = .true._LK, base = 10._RKD) )
    call disp%show("getLogSpace(0._RKD, 4._RKD, 4_IK, lopen = .true._LK, base = 10._RKD)")
    call disp%show( getLogSpace(0._RKD, 4._RKD, 4_IK, lopen = .true._LK, base = 10._RKD) )
    call disp%show("getLogSpace(0._RKD, 4._RKD, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKD)")
    call disp%show( getLogSpace(0._RKD, 4._RKD, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKD) )
    call disp%skip
    
    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace(4._RKD, 0._RKD, 4_IK)")
    call disp%show( getLogSpace(4._RKD, 0._RKD, 4_IK) )
    call disp%show("getLogSpace(4._RKD, 0._RKD, 4_IK, base = 10._RKD)")
    call disp%show( getLogSpace(4._RKD, 0._RKD, 4_IK, base = 10._RKD) )
    call disp%show("getLogSpace(4._RKD, 0._RKD, 4_IK, fopen = .true._LK, base = 10._RKD)")
    call disp%show( getLogSpace(4._RKD, 0._RKD, 4_IK, fopen = .true._LK, base = 10._RKD) )
    call disp%show("getLogSpace(4._RKD, 0._RKD, 4_IK, lopen = .true._LK, base = 10._RKD)")
    call disp%show( getLogSpace(4._RKD, 0._RKD, 4_IK, lopen = .true._LK, base = 10._RKD) )
    call disp%show("getLogSpace(4._RKD, 0._RKD, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKD)")
    call disp%show( getLogSpace(4._RKD, 0._RKD, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKD) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("128-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace(0._RKH, 4._RKH, 4_IK)")
    call disp%show( getLogSpace(0._RKH, 4._RKH, 4_IK) )
    call disp%show("getLogSpace(0._RKH, 4._RKH, 4_IK, base = 10._RKH)")
    call disp%show( getLogSpace(0._RKH, 4._RKH, 4_IK, base = 10._RKH) )
    call disp%show("getLogSpace(0._RKH, 4._RKH, 4_IK, fopen = .true._LK, base = 10._RKH)")
    call disp%show( getLogSpace(0._RKH, 4._RKH, 4_IK, fopen = .true._LK, base = 10._RKH) )
    call disp%show("getLogSpace(0._RKH, 4._RKH, 4_IK, lopen = .true._LK, base = 10._RKH)")
    call disp%show( getLogSpace(0._RKH, 4._RKH, 4_IK, lopen = .true._LK, base = 10._RKH) )
    call disp%show("getLogSpace(0._RKH, 4._RKH, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKH)")
    call disp%show( getLogSpace(0._RKH, 4._RKH, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKH) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace(4._RKH, 0._RKH, 4_IK)")
    call disp%show( getLogSpace(4._RKH, 0._RKH, 4_IK) )
    call disp%show("getLogSpace(4._RKH, 0._RKH, 4_IK, base = 10._RKH)")
    call disp%show( getLogSpace(4._RKH, 0._RKH, 4_IK, base = 10._RKH) )
    call disp%show("getLogSpace(4._RKH, 0._RKH, 4_IK, fopen = .true._LK, base = 10._RKH)")
    call disp%show( getLogSpace(4._RKH, 0._RKH, 4_IK, fopen = .true._LK, base = 10._RKH) )
    call disp%show("getLogSpace(4._RKH, 0._RKH, 4_IK, lopen = .true._LK, base = 10._RKH)")
    call disp%show( getLogSpace(4._RKH, 0._RKH, 4_IK, lopen = .true._LK, base = 10._RKH) )
    call disp%show("getLogSpace(4._RKH, 0._RKH, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKH)")
    call disp%show( getLogSpace(4._RKH, 0._RKH, 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._RKH) )
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
    call disp%show("getLogSpace((1._CKS,1._CKS), (4._CKS,4._CKS), 4_IK)")
    call disp%show( getLogSpace((1._CKS,1._CKS), (4._CKS,4._CKS), 4_IK) )
    call disp%show("getLogSpace((1._CKS,1._CKS), (4._CKS,4._CKS), 4_IK, base = 10._CKS)")
    call disp%show( getLogSpace((1._CKS,1._CKS), (4._CKS,4._CKS), 4_IK, base = 10._CKS) )
    call disp%show("getLogSpace((0._CKS,0._CKS), (4._CKS,4._CKS), 4_IK, fopen = .true._LK, base = 10._CKS)")
    call disp%show( getLogSpace((0._CKS,0._CKS), (4._CKS,4._CKS), 4_IK, fopen = .true._LK, base = 10._CKS) )
    call disp%show("getLogSpace((0._CKS,0._CKS), (4._CKS,4._CKS), 4_IK, lopen = .true._LK, base = 10._CKS)")
    call disp%show( getLogSpace((0._CKS,0._CKS), (4._CKS,4._CKS), 4_IK, lopen = .true._LK, base = 10._CKS) )
    call disp%show("getLogSpace((0._CKS,0._CKS), (4._CKS,4._CKS), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKS)")
    call disp%show( getLogSpace((0._CKS,0._CKS), (4._CKS,4._CKS), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKS) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK)")
    call disp%show( getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK) )
    call disp%show("getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK, base = 10._CKS)")
    call disp%show( getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK, base = 10._CKS) )
    call disp%show("getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK, fopen = .true._LK, base = 10._CKS)")
    call disp%show( getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK, fopen = .true._LK, base = 10._CKS) )
    call disp%show("getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK, lopen = .true._LK, base = 10._CKS)")
    call disp%show( getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK, lopen = .true._LK, base = 10._CKS) )
    call disp%show("getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKS)")
    call disp%show( getLogSpace((4._CKS,-4._CKS), (-4._CKS,+4._CKS), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKS) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("64-bit complex")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK)")
    call disp%show( getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK) )
    call disp%show("getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK, base = 10._CKD)")
    call disp%show( getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK, base = 10._CKD) )
    call disp%show("getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK, fopen = .true._LK, base = 10._CKD)")
    call disp%show( getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK, fopen = .true._LK, base = 10._CKD) )
    call disp%show("getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK, lopen = .true._LK, base = 10._CKD)")
    call disp%show( getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK, lopen = .true._LK, base = 10._CKD) )
    call disp%show("getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKD)")
    call disp%show( getLogSpace((0._CKD,0._CKD), (4._CKD,4._CKD), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKD) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK)")
    call disp%show( getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK) )
    call disp%show("getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK, base = 10._CKD)")
    call disp%show( getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK, base = 10._CKD) )
    call disp%show("getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK, fopen = .true._LK, base = 10._CKD)")
    call disp%show( getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK, fopen = .true._LK, base = 10._CKD) )
    call disp%show("getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK, lopen = .true._LK, base = 10._CKD)")
    call disp%show( getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK, lopen = .true._LK, base = 10._CKD) )
    call disp%show("getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKD)")
    call disp%show( getLogSpace((4._CKD,-4._CKD), (-4._CKD,+4._CKD), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKD) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK)")
    call disp%show( getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK) )
    call disp%show("getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK, base = 10._CKH)")
    call disp%show( getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK, base = 10._CKH) )
    call disp%show("getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK, fopen = .true._LK, base = 10._CKH)")
    call disp%show( getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK, fopen = .true._LK, base = 10._CKH) )
    call disp%show("getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK, lopen = .true._LK, base = 10._CKH)")
    call disp%show( getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK, lopen = .true._LK, base = 10._CKH) )
    call disp%show("getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKH)")
    call disp%show( getLogSpace((0._CKH,0._CKH), (4._CKH,4._CKH), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKH) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK)")
    call disp%show( getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK) )
    call disp%show("getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK, base = 10._CKH)")
    call disp%show( getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK, base = 10._CKH) )
    call disp%show("getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK, fopen = .true._LK, base = 10._CKH)")
    call disp%show( getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK, fopen = .true._LK, base = 10._CKH) )
    call disp%show("getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK, lopen = .true._LK, base = 10._CKH)")
    call disp%show( getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK, lopen = .true._LK, base = 10._CKH) )
    call disp%show("getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKH)")
    call disp%show( getLogSpace((4._CKH,-4._CKH), (-4._CKH,+4._CKH), 4_IK, fopen = .true._LK, lopen = .true._LK, base = 10._CKH) )
    call disp%skip

end program example