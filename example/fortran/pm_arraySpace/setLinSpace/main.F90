program example

    use pm_kind, only: SK
    use pm_io, only: display_type
    use pm_kind, only: IK, LK, RK32, RK64, RK128, CK32, CK64, CK128
    use pm_arraySpace, only: setLinSpace 

    implicit none

    real(RK128)     :: LogSpace_RK128(4)
    real(RK64 )     :: LogSpace_RK64 (4)
    real(RK32 )     :: LogSpace_RK32 (4)
    complex(CK128)  :: LogSpace_CK128(4)
    complex(CK64 )  :: LogSpace_CK64 (4)
    complex(CK32 )  :: LogSpace_CK32 (4)

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
    call disp%show("call setLinSpace(LogSpace_RK32, 0._RK32, 10._RK32)")
                    call setLinSpace(LogSpace_RK32, 0._RK32, 10._RK32)
    call disp%show("LogSpace_RK32")
    call disp%show( LogSpace_RK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK32, 0._RK32, 10._RK32, fopen = .true._LK)")
                    call setLinSpace(LogSpace_RK32, 0._RK32, 10._RK32, fopen = .true._LK)
    call disp%show("LogSpace_RK32")
    call disp%show( LogSpace_RK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK32, 0._RK32, 10._RK32, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK32, 0._RK32, 10._RK32, lopen = .true._LK)
    call disp%show("LogSpace_RK32")
    call disp%show( LogSpace_RK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK32, 0._RK32, 10._RK32, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK32, 0._RK32, 10._RK32, fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_RK32")
    call disp%show( LogSpace_RK32 )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK32, 10._RK32, 0._RK32)")
                    call setLinSpace(LogSpace_RK32, 10._RK32, 0._RK32)
    call disp%show("LogSpace_RK32")
    call disp%show( LogSpace_RK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK32, 10._RK32, 0._RK32, fopen = .true._LK)")
                    call setLinSpace(LogSpace_RK32, 10._RK32, 0._RK32, fopen = .true._LK)
    call disp%show("LogSpace_RK32")
    call disp%show( LogSpace_RK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK32, 10._RK32, 0._RK32, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK32, 10._RK32, 0._RK32, lopen = .true._LK)
    call disp%show("LogSpace_RK32")
    call disp%show( LogSpace_RK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK32, 10._RK32, 0._RK32, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK32, 10._RK32, 0._RK32, fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_RK32")
    call disp%show( LogSpace_RK32 )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!64-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK64, 0._RK64, 10._RK64)")
                    call setLinSpace(LogSpace_RK64, 0._RK64, 10._RK64)
    call disp%show("LogSpace_RK64")
    call disp%show( LogSpace_RK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK64, 0._RK64, 10._RK64, fopen = .true._LK)")
                    call setLinSpace(LogSpace_RK64, 0._RK64, 10._RK64, fopen = .true._LK)
    call disp%show("LogSpace_RK64")
    call disp%show( LogSpace_RK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK64, 0._RK64, 10._RK64, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK64, 0._RK64, 10._RK64, lopen = .true._LK)
    call disp%show("LogSpace_RK64")
    call disp%show( LogSpace_RK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK64, 0._RK64, 10._RK64, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK64, 0._RK64, 10._RK64, fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_RK64")
    call disp%show( LogSpace_RK64 )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK64, 10._RK64, 0._RK64)")
                    call setLinSpace(LogSpace_RK64, 10._RK64, 0._RK64)
    call disp%show("LogSpace_RK64")
    call disp%show( LogSpace_RK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK64, 10._RK64, 0._RK64, fopen = .true._LK)")
                    call setLinSpace(LogSpace_RK64, 10._RK64, 0._RK64, fopen = .true._LK)
    call disp%show("LogSpace_RK64")
    call disp%show( LogSpace_RK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK64, 10._RK64, 0._RK64, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK64, 10._RK64, 0._RK64, lopen = .true._LK)
    call disp%show("LogSpace_RK64")
    call disp%show( LogSpace_RK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK64, 10._RK64, 0._RK64, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK64, 10._RK64, 0._RK64, fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_RK64")
    call disp%show( LogSpace_RK64 )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!128-bit real")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK128, 0._RK128, 10._RK128)")
                    call setLinSpace(LogSpace_RK128, 0._RK128, 10._RK128)
    call disp%show("LogSpace_RK128")
    call disp%show( LogSpace_RK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK128, 0._RK128, 10._RK128, fopen = .true._LK)")
                    call setLinSpace(LogSpace_RK128, 0._RK128, 10._RK128, fopen = .true._LK)
    call disp%show("LogSpace_RK128")
    call disp%show( LogSpace_RK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK128, 0._RK128, 10._RK128, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK128, 0._RK128, 10._RK128, lopen = .true._LK)
    call disp%show("LogSpace_RK128")
    call disp%show( LogSpace_RK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK128, 0._RK128, 10._RK128, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK128, 0._RK128, 10._RK128, fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_RK128")
    call disp%show( LogSpace_RK128 )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK128, 10._RK128, 0._RK128)")
                    call setLinSpace(LogSpace_RK128, 10._RK128, 0._RK128)
    call disp%show("LogSpace_RK128")
    call disp%show( LogSpace_RK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK128, 10._RK128, 0._RK128, fopen = .true._LK)")
                    call setLinSpace(LogSpace_RK128, 10._RK128, 0._RK128, fopen = .true._LK)
    call disp%show("LogSpace_RK128")
    call disp%show( LogSpace_RK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK128, 10._RK128, 0._RK128, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK128, 10._RK128, 0._RK128, lopen = .true._LK)
    call disp%show("LogSpace_RK128")
    call disp%show( LogSpace_RK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_RK128, 10._RK128, 0._RK128, fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_RK128, 10._RK128, 0._RK128, fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_RK128")
    call disp%show( LogSpace_RK128 )
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
    call disp%show("call setLinSpace(LogSpace_CK32, (0._CK32,0._CK32), (10._CK32,10._CK32))")
                    call setLinSpace(LogSpace_CK32, (0._CK32,0._CK32), (10._CK32,10._CK32))
    call disp%show("LogSpace_CK32")
    call disp%show( LogSpace_CK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK32, (0._CK32,0._CK32), (10._CK32,10._CK32), fopen = .true._LK)")
                    call setLinSpace(LogSpace_CK32, (0._CK32,0._CK32), (10._CK32,10._CK32), fopen = .true._LK)
    call disp%show("LogSpace_CK32")
    call disp%show( LogSpace_CK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK32, (0._CK32,0._CK32), (10._CK32,10._CK32), lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK32, (0._CK32,0._CK32), (10._CK32,10._CK32), lopen = .true._LK)
    call disp%show("LogSpace_CK32")
    call disp%show( LogSpace_CK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK32, (0._CK32,0._CK32), (10._CK32,10._CK32), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK32, (0._CK32,0._CK32), (10._CK32,10._CK32), fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_CK32")
    call disp%show( LogSpace_CK32 )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK32, (10._CK32,-10._CK32), (0._CK32,0._CK32))")
                    call setLinSpace(LogSpace_CK32, (10._CK32,-10._CK32), (0._CK32,0._CK32))
    call disp%show("LogSpace_CK32")
    call disp%show( LogSpace_CK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK32, (10._CK32,-10._CK32), (0._CK32,0._CK32), fopen = .true._LK)")
                    call setLinSpace(LogSpace_CK32, (10._CK32,-10._CK32), (0._CK32,0._CK32), fopen = .true._LK)
    call disp%show("LogSpace_CK32")
    call disp%show( LogSpace_CK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK32, (10._CK32,-10._CK32), (0._CK32,0._CK32), lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK32, (10._CK32,-10._CK32), (0._CK32,0._CK32), lopen = .true._LK)
    call disp%show("LogSpace_CK32")
    call disp%show( LogSpace_CK32 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK32, (10._CK32,-10._CK32), (0._CK32,0._CK32), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK32, (10._CK32,-10._CK32), (0._CK32,0._CK32), fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_CK32")
    call disp%show( LogSpace_CK32 )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK64, (0._CK64,0._CK64), (10._CK64,10._CK64))")
                    call setLinSpace(LogSpace_CK64, (0._CK64,0._CK64), (10._CK64,10._CK64))
    call disp%show("LogSpace_CK64")
    call disp%show( LogSpace_CK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK64, (0._CK64,0._CK64), (10._CK64,10._CK64), fopen = .true._LK)")
                    call setLinSpace(LogSpace_CK64, (0._CK64,0._CK64), (10._CK64,10._CK64), fopen = .true._LK)
    call disp%show("LogSpace_CK64")
    call disp%show( LogSpace_CK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK64, (0._CK64,0._CK64), (10._CK64,10._CK64), lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK64, (0._CK64,0._CK64), (10._CK64,10._CK64), lopen = .true._LK)
    call disp%show("LogSpace_CK64")
    call disp%show( LogSpace_CK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK64, (0._CK64,0._CK64), (10._CK64,10._CK64), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK64, (0._CK64,0._CK64), (10._CK64,10._CK64), fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_CK64")
    call disp%show( LogSpace_CK64 )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK64, (10._CK64,-10._CK64), (0._CK64,0._CK64))")
                    call setLinSpace(LogSpace_CK64, (10._CK64,-10._CK64), (0._CK64,0._CK64))
    call disp%show("LogSpace_CK64")
    call disp%show( LogSpace_CK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK64, (10._CK64,-10._CK64), (0._CK64,0._CK64), fopen = .true._LK)")
                    call setLinSpace(LogSpace_CK64, (10._CK64,-10._CK64), (0._CK64,0._CK64), fopen = .true._LK)
    call disp%show("LogSpace_CK64")
    call disp%show( LogSpace_CK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK64, (10._CK64,-10._CK64), (0._CK64,0._CK64), lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK64, (10._CK64,-10._CK64), (0._CK64,0._CK64), lopen = .true._LK)
    call disp%show("LogSpace_CK64")
    call disp%show( LogSpace_CK64 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK64, (10._CK64,-10._CK64), (0._CK64,0._CK64), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK64, (10._CK64,-10._CK64), (0._CK64,0._CK64), fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_CK64")
    call disp%show( LogSpace_CK64 )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK128, (0._CK128,0._CK128), (10._CK128,10._CK128))")
                    call setLinSpace(LogSpace_CK128, (0._CK128,0._CK128), (10._CK128,10._CK128))
    call disp%show("LogSpace_CK128")
    call disp%show( LogSpace_CK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK128, (0._CK128,0._CK128), (10._CK128,10._CK128), fopen = .true._LK)")
                    call setLinSpace(LogSpace_CK128, (0._CK128,0._CK128), (10._CK128,10._CK128), fopen = .true._LK)
    call disp%show("LogSpace_CK128")
    call disp%show( LogSpace_CK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK128, (0._CK128,0._CK128), (10._CK128,10._CK128), lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK128, (0._CK128,0._CK128), (10._CK128,10._CK128), lopen = .true._LK)
    call disp%show("LogSpace_CK128")
    call disp%show( LogSpace_CK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK128, (0._CK128,0._CK128), (10._CK128,10._CK128), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK128, (0._CK128,0._CK128), (10._CK128,10._CK128), fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_CK128")
    call disp%show( LogSpace_CK128 )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK128, (10._CK128,-10._CK128), (0._CK128,0._CK128))")
                    call setLinSpace(LogSpace_CK128, (10._CK128,-10._CK128), (0._CK128,0._CK128))
    call disp%show("LogSpace_CK128")
    call disp%show( LogSpace_CK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK128, (10._CK128,-10._CK128), (0._CK128,0._CK128), fopen = .true._LK)")
                    call setLinSpace(LogSpace_CK128, (10._CK128,-10._CK128), (0._CK128,0._CK128), fopen = .true._LK)
    call disp%show("LogSpace_CK128")
    call disp%show( LogSpace_CK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK128, (10._CK128,-10._CK128), (0._CK128,0._CK128), lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK128, (10._CK128,-10._CK128), (0._CK128,0._CK128), lopen = .true._LK)
    call disp%show("LogSpace_CK128")
    call disp%show( LogSpace_CK128 )
    call disp%skip
    call disp%show("call setLinSpace(LogSpace_CK128, (10._CK128,-10._CK128), (0._CK128,0._CK128), fopen = .true._LK, lopen = .true._LK)")
                    call setLinSpace(LogSpace_CK128, (10._CK128,-10._CK128), (0._CK128,0._CK128), fopen = .true._LK, lopen = .true._LK)
    call disp%show("LogSpace_CK128")
    call disp%show( LogSpace_CK128 )
    call disp%skip

end program example