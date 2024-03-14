program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RK32, RK64, RK128 ! all processor types and kinds are supported.
    use pm_kind, only: CK32, CK64, CK128 ! all processor types and kinds are supported.
    use pm_except, only: isInfNeg, setInfNeg
    use pm_io, only: display_type

    implicit none

    real(RK128)     :: X_RK128(3)
    real(RK64 )     :: X_RK64 (3)
    real(RK32 )     :: X_RK32 (3)
    complex(CK128)  :: X_CK128(3)
    complex(CK64 )  :: X_CK64 (3)
    complex(CK32 )  :: X_CK32 (3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setInfNeg(X_RK128)
    call setInfNeg(X_RK64 )
    call setInfNeg(X_RK32 )

    call setInfNeg(X_CK128)
    call setInfNeg(X_CK64 )
    call setInfNeg(X_CK32 )

    X_RK128(2) = 0._RK128
    X_RK64 (2) = 0._RK64
    X_RK32 (2) = 0._RK32

    X_CK128(2) = (0._CK128, 0._CK128)
    X_CK64 (2) = (0._CK64 , 0._CK64 )
    X_CK32 (2) = (0._CK32 , 0._CK32 )

    X_CK128(2) = cmplx(0._CK128     , X_RK128(1), CK128)
    X_CK64 (2) = cmplx(-X_RK64(1)   , 0._CK64   , CK64 )
    X_CK32 (2) = cmplx(0._CK32      , 0._CK32   , CK32 )

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Generate real IEEE-compliant negative infinity.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!32-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_RK32(1)")
    call disp%show( X_RK32(1) )
    call disp%show("isInfNeg(X_RK32(1))")
    call disp%show( isInfNeg(X_RK32(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_RK32")
    call disp%show( X_RK32 )
    call disp%show("isInfNeg(X_RK32)")
    call disp%show( isInfNeg(X_RK32) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!64-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_RK64(1)")
    call disp%show( X_RK64(1) )
    call disp%show("isInfNeg(X_RK64(1))")
    call disp%show( isInfNeg(X_RK64(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_RK64")
    call disp%show( X_RK64 )
    call disp%show("isInfNeg(X_RK64)")
    call disp%show( isInfNeg(X_RK64) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!128-bit real")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_RK128(1)")
    call disp%show( X_RK128(1) )
    call disp%show("isInfNeg(X_RK128(1))")
    call disp%show( isInfNeg(X_RK128(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_RK128")
    call disp%show( X_RK128 )
    call disp%show("isInfNeg(X_RK128)")
    call disp%show( isInfNeg(X_RK128) )
    call disp%skip


    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Generate complex IEEE-compliant negative infinity.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip


    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!32-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_CK32(1)")
    call disp%show( X_CK32(1) )
    call disp%show("isInfNeg(X_CK32(1))")
    call disp%show( isInfNeg(X_CK32(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_CK32")
    call disp%show( X_CK32 )
    call disp%show("isInfNeg(X_CK32)")
    call disp%show( isInfNeg(X_CK32) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_CK64(1)")
    call disp%show( X_CK64(1) )
    call disp%show("isInfNeg(X_CK64(1))")
    call disp%show( isInfNeg(X_CK64(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_CK64")
    call disp%show( X_CK64 )
    call disp%show("isInfNeg(X_CK64)")
    call disp%show( isInfNeg(X_CK64) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_CK128(1)")
    call disp%show( X_CK128(1) )
    call disp%show("isInfNeg(X_CK128(1))")
    call disp%show( isInfNeg(X_CK128(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_CK128")
    call disp%show( X_CK128 )
    call disp%show("isInfNeg(X_CK128)")
    call disp%show( isInfNeg(X_CK128) )
    call disp%skip

end program example