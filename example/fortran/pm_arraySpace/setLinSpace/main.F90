program example

    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_kind, only: IK, LK, RKS, RKD, RKH, CKS, CKD, CKH
    use pm_arraySpace, only: setLinSpace 

    implicit none

    real(RKH)       :: logSpace_RKH(4)
    real(RKD)       :: logSpace_RKD(4)
    real(RKS)       :: logSpace_RKS(4)
    complex(CKH)    :: logSpace_CKH(4)
    complex(CKD)    :: logSpace_CKD(4)
    complex(CKS)    :: logSpace_CKS(4)

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
    call disp%show("call setLinSpace(logSpace_RKS, 0._RKS, 10._RKS)")
                    call setLinSpace(logSpace_RKS, 0._RKS, 10._RKS)
    call disp%show("logSpace_RKS")
    call disp%show( logSpace_RKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKS, 0._RKS, 10._RKS, fopen = .true._LK)")
                    call setLinSpace(logSpace_RKS, 0._RKS, 10._RKS, fopen = .true._LK)
    call disp%show("logSpace_RKS")
    call disp%show( logSpace_RKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKS, 0._RKS, 10._RKS, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKS, 0._RKS, 10._RKS, lopen = .true._LK)
    call disp%show("logSpace_RKS")
    call disp%show( logSpace_RKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKS, 0._RKS, 10._RKS, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKS, 0._RKS, 10._RKS, fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_RKS")
    call disp%show( logSpace_RKS )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKS, 10._RKS, 0._RKS)")
                    call setLinSpace(logSpace_RKS, 10._RKS, 0._RKS)
    call disp%show("logSpace_RKS")
    call disp%show( logSpace_RKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKS, 10._RKS, 0._RKS, fopen = .true._LK)")
                    call setLinSpace(logSpace_RKS, 10._RKS, 0._RKS, fopen = .true._LK)
    call disp%show("logSpace_RKS")
    call disp%show( logSpace_RKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKS, 10._RKS, 0._RKS, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKS, 10._RKS, 0._RKS, lopen = .true._LK)
    call disp%show("logSpace_RKS")
    call disp%show( logSpace_RKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKS, 10._RKS, 0._RKS, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKS, 10._RKS, 0._RKS, fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_RKS")
    call disp%show( logSpace_RKS )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!64-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKD, 0._RKD, 10._RKD)")
                    call setLinSpace(logSpace_RKD, 0._RKD, 10._RKD)
    call disp%show("logSpace_RKD")
    call disp%show( logSpace_RKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKD, 0._RKD, 10._RKD, fopen = .true._LK)")
                    call setLinSpace(logSpace_RKD, 0._RKD, 10._RKD, fopen = .true._LK)
    call disp%show("logSpace_RKD")
    call disp%show( logSpace_RKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKD, 0._RKD, 10._RKD, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKD, 0._RKD, 10._RKD, lopen = .true._LK)
    call disp%show("logSpace_RKD")
    call disp%show( logSpace_RKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKD, 0._RKD, 10._RKD, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKD, 0._RKD, 10._RKD, fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_RKD")
    call disp%show( logSpace_RKD )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKD, 10._RKD, 0._RKD)")
                    call setLinSpace(logSpace_RKD, 10._RKD, 0._RKD)
    call disp%show("logSpace_RKD")
    call disp%show( logSpace_RKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKD, 10._RKD, 0._RKD, fopen = .true._LK)")
                    call setLinSpace(logSpace_RKD, 10._RKD, 0._RKD, fopen = .true._LK)
    call disp%show("logSpace_RKD")
    call disp%show( logSpace_RKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKD, 10._RKD, 0._RKD, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKD, 10._RKD, 0._RKD, lopen = .true._LK)
    call disp%show("logSpace_RKD")
    call disp%show( logSpace_RKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKD, 10._RKD, 0._RKD, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKD, 10._RKD, 0._RKD, fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_RKD")
    call disp%show( logSpace_RKD )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!128-bit real")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKH, 0._RKH, 10._RKH)")
                    call setLinSpace(logSpace_RKH, 0._RKH, 10._RKH)
    call disp%show("logSpace_RKH")
    call disp%show( logSpace_RKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKH, 0._RKH, 10._RKH, fopen = .true._LK)")
                    call setLinSpace(logSpace_RKH, 0._RKH, 10._RKH, fopen = .true._LK)
    call disp%show("logSpace_RKH")
    call disp%show( logSpace_RKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKH, 0._RKH, 10._RKH, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKH, 0._RKH, 10._RKH, lopen = .true._LK)
    call disp%show("logSpace_RKH")
    call disp%show( logSpace_RKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKH, 0._RKH, 10._RKH, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKH, 0._RKH, 10._RKH, fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_RKH")
    call disp%show( logSpace_RKH )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKH, 10._RKH, 0._RKH)")
                    call setLinSpace(logSpace_RKH, 10._RKH, 0._RKH)
    call disp%show("logSpace_RKH")
    call disp%show( logSpace_RKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKH, 10._RKH, 0._RKH, fopen = .true._LK)")
                    call setLinSpace(logSpace_RKH, 10._RKH, 0._RKH, fopen = .true._LK)
    call disp%show("logSpace_RKH")
    call disp%show( logSpace_RKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKH, 10._RKH, 0._RKH, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKH, 10._RKH, 0._RKH, lopen = .true._LK)
    call disp%show("logSpace_RKH")
    call disp%show( logSpace_RKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_RKH, 10._RKH, 0._RKH, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_RKH, 10._RKH, 0._RKH, fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_RKH")
    call disp%show( logSpace_RKH )
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
    call disp%show("call setLinSpace(logSpace_CKS, (0._CKS,0._CKS), (10._CKS,10._CKS))")
                    call setLinSpace(logSpace_CKS, (0._CKS,0._CKS), (10._CKS,10._CKS))
    call disp%show("logSpace_CKS")
    call disp%show( logSpace_CKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKS, (0._CKS,0._CKS), (10._CKS,10._CKS), fopen = .true._LK)")
                    call setLinSpace(logSpace_CKS, (0._CKS,0._CKS), (10._CKS,10._CKS), fopen = .true._LK)
    call disp%show("logSpace_CKS")
    call disp%show( logSpace_CKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKS, (0._CKS,0._CKS), (10._CKS,10._CKS), lopen = .true._LK)")
                    call setLinSpace(logSpace_CKS, (0._CKS,0._CKS), (10._CKS,10._CKS), lopen = .true._LK)
    call disp%show("logSpace_CKS")
    call disp%show( logSpace_CKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKS, (0._CKS,0._CKS), (10._CKS,10._CKS), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_CKS, (0._CKS,0._CKS), (10._CKS,10._CKS), fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_CKS")
    call disp%show( logSpace_CKS )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKS, (10._CKS,-10._CKS), (0._CKS,0._CKS))")
                    call setLinSpace(logSpace_CKS, (10._CKS,-10._CKS), (0._CKS,0._CKS))
    call disp%show("logSpace_CKS")
    call disp%show( logSpace_CKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKS, (10._CKS,-10._CKS), (0._CKS,0._CKS), fopen = .true._LK)")
                    call setLinSpace(logSpace_CKS, (10._CKS,-10._CKS), (0._CKS,0._CKS), fopen = .true._LK)
    call disp%show("logSpace_CKS")
    call disp%show( logSpace_CKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKS, (10._CKS,-10._CKS), (0._CKS,0._CKS), lopen = .true._LK)")
                    call setLinSpace(logSpace_CKS, (10._CKS,-10._CKS), (0._CKS,0._CKS), lopen = .true._LK)
    call disp%show("logSpace_CKS")
    call disp%show( logSpace_CKS )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKS, (10._CKS,-10._CKS), (0._CKS,0._CKS), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_CKS, (10._CKS,-10._CKS), (0._CKS,0._CKS), fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_CKS")
    call disp%show( logSpace_CKS )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKD, (0._CKD,0._CKD), (10._CKD,10._CKD))")
                    call setLinSpace(logSpace_CKD, (0._CKD,0._CKD), (10._CKD,10._CKD))
    call disp%show("logSpace_CKD")
    call disp%show( logSpace_CKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKD, (0._CKD,0._CKD), (10._CKD,10._CKD), fopen = .true._LK)")
                    call setLinSpace(logSpace_CKD, (0._CKD,0._CKD), (10._CKD,10._CKD), fopen = .true._LK)
    call disp%show("logSpace_CKD")
    call disp%show( logSpace_CKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKD, (0._CKD,0._CKD), (10._CKD,10._CKD), lopen = .true._LK)")
                    call setLinSpace(logSpace_CKD, (0._CKD,0._CKD), (10._CKD,10._CKD), lopen = .true._LK)
    call disp%show("logSpace_CKD")
    call disp%show( logSpace_CKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKD, (0._CKD,0._CKD), (10._CKD,10._CKD), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_CKD, (0._CKD,0._CKD), (10._CKD,10._CKD), fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_CKD")
    call disp%show( logSpace_CKD )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKD, (10._CKD,-10._CKD), (0._CKD,0._CKD))")
                    call setLinSpace(logSpace_CKD, (10._CKD,-10._CKD), (0._CKD,0._CKD))
    call disp%show("logSpace_CKD")
    call disp%show( logSpace_CKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKD, (10._CKD,-10._CKD), (0._CKD,0._CKD), fopen = .true._LK)")
                    call setLinSpace(logSpace_CKD, (10._CKD,-10._CKD), (0._CKD,0._CKD), fopen = .true._LK)
    call disp%show("logSpace_CKD")
    call disp%show( logSpace_CKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKD, (10._CKD,-10._CKD), (0._CKD,0._CKD), lopen = .true._LK)")
                    call setLinSpace(logSpace_CKD, (10._CKD,-10._CKD), (0._CKD,0._CKD), lopen = .true._LK)
    call disp%show("logSpace_CKD")
    call disp%show( logSpace_CKD )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKD, (10._CKD,-10._CKD), (0._CKD,0._CKD), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_CKD, (10._CKD,-10._CKD), (0._CKD,0._CKD), fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_CKD")
    call disp%show( logSpace_CKD )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKH, (0._CKH,0._CKH), (10._CKH,10._CKH))")
                    call setLinSpace(logSpace_CKH, (0._CKH,0._CKH), (10._CKH,10._CKH))
    call disp%show("logSpace_CKH")
    call disp%show( logSpace_CKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKH, (0._CKH,0._CKH), (10._CKH,10._CKH), fopen = .true._LK)")
                    call setLinSpace(logSpace_CKH, (0._CKH,0._CKH), (10._CKH,10._CKH), fopen = .true._LK)
    call disp%show("logSpace_CKH")
    call disp%show( logSpace_CKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKH, (0._CKH,0._CKH), (10._CKH,10._CKH), lopen = .true._LK)")
                    call setLinSpace(logSpace_CKH, (0._CKH,0._CKH), (10._CKH,10._CKH), lopen = .true._LK)
    call disp%show("logSpace_CKH")
    call disp%show( logSpace_CKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKH, (0._CKH,0._CKH), (10._CKH,10._CKH), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_CKH, (0._CKH,0._CKH), (10._CKH,10._CKH), fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_CKH")
    call disp%show( logSpace_CKH )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKH, (10._CKH,-10._CKH), (0._CKH,0._CKH))")
                    call setLinSpace(logSpace_CKH, (10._CKH,-10._CKH), (0._CKH,0._CKH))
    call disp%show("logSpace_CKH")
    call disp%show( logSpace_CKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKH, (10._CKH,-10._CKH), (0._CKH,0._CKH), fopen = .true._LK)")
                    call setLinSpace(logSpace_CKH, (10._CKH,-10._CKH), (0._CKH,0._CKH), fopen = .true._LK)
    call disp%show("logSpace_CKH")
    call disp%show( logSpace_CKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKH, (10._CKH,-10._CKH), (0._CKH,0._CKH), lopen = .true._LK)")
                    call setLinSpace(logSpace_CKH, (10._CKH,-10._CKH), (0._CKH,0._CKH), lopen = .true._LK)
    call disp%show("logSpace_CKH")
    call disp%show( logSpace_CKH )
    call disp%skip
    call disp%show("call setLinSpace(logSpace_CKH, (10._CKH,-10._CKH), (0._CKH,0._CKH), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(logSpace_CKH, (10._CKH,-10._CKH), (0._CKH,0._CKH), fopen = .true._LK, lopen = .true._LK)
    call disp%show("logSpace_CKH")
    call disp%show( logSpace_CKH )
    call disp%skip

end program example