program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKS, RKD, RKH ! all processor types and kinds are supported.
    use pm_kind, only: CKS, CKD, CKH ! all processor types and kinds are supported.
    use pm_except, only: isInfNeg, setInfNeg
    use pm_io, only: display_type

    implicit none

    real(RKH)     :: X_RKH(3)
    real(RKD)     :: X_RKD(3)
    real(RKS)     :: X_RKS(3)
    complex(CKH)  :: X_CKH(3)
    complex(CKD)  :: X_CKD(3)
    complex(CKS)  :: X_CKS(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setInfNeg(X_RKH)
    call setInfNeg(X_RKD)
    call setInfNeg(X_RKS)

    call setInfNeg(X_CKH)
    call setInfNeg(X_CKD)
    call setInfNeg(X_CKS)

    X_RKH(2) = 0._RKH
    X_RKD(2) = 0._RKD
    X_RKS(2) = 0._RKS

    X_CKH(2) = (0._CKH, 0._CKH)
    X_CKD(2) = (0._CKD, 0._CKD)
    X_CKS(2) = (0._CKS, 0._CKS)

    X_CKH(2) = cmplx(0._CKH, X_RKH(1), CKH)
    X_CKD(2) = cmplx(-X_RKD(1), 0._CKD, CKD)
    X_CKS(2) = cmplx(0._CKS, 0._CKS, CKS)

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
    call disp%show("X_RKS(1)")
    call disp%show( X_RKS(1) )
    call disp%show("isInfNeg(X_RKS(1))")
    call disp%show( isInfNeg(X_RKS(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_RKS")
    call disp%show( X_RKS )
    call disp%show("isInfNeg(X_RKS)")
    call disp%show( isInfNeg(X_RKS) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!64-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_RKD(1)")
    call disp%show( X_RKD(1) )
    call disp%show("isInfNeg(X_RKD(1))")
    call disp%show( isInfNeg(X_RKD(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_RKD")
    call disp%show( X_RKD)
    call disp%show("isInfNeg(X_RKD)")
    call disp%show( isInfNeg(X_RKD) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!128-bit real")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_RKH(1)")
    call disp%show( X_RKH(1) )
    call disp%show("isInfNeg(X_RKH(1))")
    call disp%show( isInfNeg(X_RKH(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_RKH")
    call disp%show( X_RKH )
    call disp%show("isInfNeg(X_RKH)")
    call disp%show( isInfNeg(X_RKH) )
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
    call disp%show("X_CKS(1)")
    call disp%show( X_CKS(1) )
    call disp%show("isInfNeg(X_CKS(1))")
    call disp%show( isInfNeg(X_CKS(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_CKS")
    call disp%show( X_CKS )
    call disp%show("isInfNeg(X_CKS)")
    call disp%show( isInfNeg(X_CKS) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_CKD(1)")
    call disp%show( X_CKD(1) )
    call disp%show("isInfNeg(X_CKD(1))")
    call disp%show( isInfNeg(X_CKD(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_CKD")
    call disp%show( X_CKD )
    call disp%show("isInfNeg(X_CKD)")
    call disp%show( isInfNeg(X_CKD) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_CKH(1)")
    call disp%show( X_CKH(1) )
    call disp%show("isInfNeg(X_CKH(1))")
    call disp%show( isInfNeg(X_CKH(1)) )
    call disp%skip

    call disp%skip
    call disp%show("X_CKH")
    call disp%show( X_CKH )
    call disp%show("isInfNeg(X_CKH)")
    call disp%show( isInfNeg(X_CKH) )
    call disp%skip

end program example