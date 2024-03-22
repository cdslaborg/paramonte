program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKS, RKD, RKH ! all processor types and kinds are supported.
    use pm_kind, only: CKS, CKD, CKH ! all processor types and kinds are supported.
    use pm_io, only: display_type
    use pm_except, only: getInfPos

    implicit none

    real(RKH)     :: X_RKH(3)
    real(RKD)     :: X_RKD(3)
    real(RKS)     :: X_RKS(3)
    complex(CKH)  :: X_CKH(3)
    complex(CKD)  :: X_CKD(3)
    complex(CKS)  :: X_CKS(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Generate real positive Infinity.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!32-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_RKS(1) = getInfPos(mold = X_RKS(1))")
                    X_RKS(1) = getInfPos(mold = X_RKS(1))
    call disp%show("X_RKS(1)")
    call disp%show( X_RKS(1) )
    call disp%skip

    call disp%skip
    call disp%show("X_RKS = getInfPos(mold = X_RKS)")
                    X_RKS = getInfPos(mold = X_RKS)
    call disp%show("X_RKS")
    call disp%show( X_RKS )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%")
    call disp%show("!64-bit real")
    call disp%show("!%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_RKD(1) = getInfPos(mold = X_RKD(1))")
                    X_RKD(1) = getInfPos(mold = X_RKD(1))
    call disp%show("X_RKD(1)")
    call disp%show( X_RKD(1) )
    call disp%skip

    call disp%skip
    call disp%show("X_RKD = getInfPos(mold = X_RKD)")
                    X_RKD = getInfPos(mold = X_RKD)
    call disp%show("X_RKD")
    call disp%show( X_RKD )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%")
    call disp%show("!128-bit real")
    call disp%show("!%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_RKH(1) = getInfPos(mold = X_RKH(1))")
                    X_RKH(1) = getInfPos(mold = X_RKH(1))
    call disp%show("X_RKH(1)")
    call disp%show( X_RKH(1) )
    call disp%skip

    call disp%skip
    call disp%show("X_RKH = getInfPos(mold = X_RKH)")
                    X_RKH = getInfPos(mold = X_RKH)
    call disp%show("X_RKH")
    call disp%show( X_RKH )
    call disp%skip


    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Generate complex positive Infinity.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip


    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!32-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_CKS(1) = getInfPos(mold = X_CKS(1))")
                    X_CKS(1) = getInfPos(mold = X_CKS(1))
    call disp%show("X_CKS(1)")
    call disp%show( X_CKS(1) )
    call disp%skip

    call disp%skip
    call disp%show("X_CKS = getInfPos(mold = X_CKS)")
                    X_CKS = getInfPos(mold = X_CKS)
    call disp%show("X_CKS")
    call disp%show( X_CKS )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit complex")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_CKD(1) = getInfPos(mold = X_CKD(1))")
                    X_CKD(1) = getInfPos(mold = X_CKD(1))
    call disp%show("X_CKD(1)")
    call disp%show( X_CKD(1) )
    call disp%skip

    call disp%skip
    call disp%show("X_CKD = getInfPos(mold = X_CKD)")
                    X_CKD = getInfPos(mold = X_CKD)
    call disp%show("X_CKD")
    call disp%show( X_CKD )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!128-bit complex")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("X_CKH(1) = getInfPos(mold = X_CKH(1))")
                    X_CKH(1) = getInfPos(mold = X_CKH(1))
    call disp%show("X_CKH(1)")
    call disp%show( X_CKH(1) )
    call disp%skip

    call disp%skip
    call disp%show("X_CKH = getInfPos(mold = X_CKH)")
                    X_CKH = getInfPos(mold = X_CKH)
    call disp%show("X_CKH")
    call disp%show( X_CKH )
    call disp%skip

end program example