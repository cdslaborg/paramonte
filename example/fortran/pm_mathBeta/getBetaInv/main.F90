program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_mathBeta, only: getBetaInv
    use pm_distUnif, only: getUnifRand
    use pm_mathBeta, only: getBetaInc

    implicit none

    real(RKC), allocatable :: alpha, beta, betaInc(:), betaInv(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("betaInc = [0.5_RKC]")
                    betaInc = [0.5_RKC]
    call disp%show("alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)")
                    alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)
    call disp%show("[alpha, beta]")
    call disp%show( [alpha, beta] )
    call disp%show("betaInv = getBetaInv(betaInc, alpha, beta)")
                    betaInv = getBetaInv(betaInc, alpha, beta)
    call disp%show("betaInv")
    call disp%show( betaInv )
    call disp%show("getBetaInc(betaInv, alpha, beta)")
    call disp%show( getBetaInc(betaInv, alpha, beta) )
    call disp%skip()

    call disp%skip()
    call disp%show("betaInc = [0.5_RKC]")
                    betaInc = [0.5_RKC]
    call disp%show("alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)")
                    alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)
    call disp%show("[alpha, beta]")
    call disp%show( [alpha, beta] )
    call disp%show("betaInv = getBetaInv(betaInc, alpha, beta)")
                    betaInv = getBetaInv(betaInc, alpha, beta)
    call disp%show("betaInv")
    call disp%show( betaInv )
    call disp%show("getBetaInc(betaInv, alpha, beta)")
    call disp%show( getBetaInc(betaInv, alpha, beta) )
    call disp%skip()

    call disp%skip()
    call disp%show("betaInc = [1._RKC - epsilon(1._RKC)]")
                    betaInc = [1._RKC - epsilon(1._RKC)]
    call disp%show("alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)")
                    alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)
    call disp%show("[alpha, beta]")
    call disp%show( [alpha, beta] )
    call disp%show("betaInv = getBetaInv(betaInc, alpha, beta)")
                    betaInv = getBetaInv(betaInc, alpha, beta)
    call disp%show("betaInv")
    call disp%show( betaInv )
    call disp%show("getBetaInc(betaInv, alpha, beta)")
    call disp%show( getBetaInc(betaInv, alpha, beta) )
    call disp%skip()

    call disp%skip()
    call disp%show("betaInc = [0._RKC, .5_RKC, 1._RKC]")
                    betaInc = [0._RKC, .5_RKC, 1._RKC]
    call disp%show("alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)")
                    alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)
    call disp%show("[alpha, beta]")
    call disp%show( [alpha, beta] )
    call disp%show("betaInv = getBetaInv(betaInc, alpha, beta)")
                    betaInv = getBetaInv(betaInc, alpha, beta)
    call disp%show("betaInv")
    call disp%show( betaInv )
    call disp%show("getBetaInc(betaInv, alpha, beta)")
    call disp%show( getBetaInc(betaInv, alpha, beta) )
    call disp%skip()

    call disp%skip()
    call disp%show("betaInc = [0._RKC, .5_RKC, 1._RKC]")
                    betaInc = [0._RKC, .5_RKC, 1._RKC]
    call disp%show("betaInv = getBetaInv(betaInc, alpha = [.1_RKC, 1._RKC, 10._RKC], beta = 3._RKC)")
                    betaInv = getBetaInv(betaInc, alpha = [.1_RKC, 1._RKC, 10._RKC], beta = 3._RKC)
    call disp%show("betaInv")
    call disp%show( betaInv )
    call disp%show("getBetaInc(betaInv, alpha = [.1_RKC, 1._RKC, 10._RKC], beta = 3._RKC)")
    call disp%show( getBetaInc(betaInv, alpha = [.1_RKC, 1._RKC, 10._RKC], beta = 3._RKC) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_kind, only: RKC => RKH, RKH
        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 1000
        real(RKC)   , parameter :: alpha(*) = [0.1_RKC, 10._RKC, 1._RKC, 0.1_RKC, 10._RKC], beta(*) = [0.1_RKC, 0.1_RKC, 1._RKC, 10._RKC, 10._RKC]
        real(RKC)   :: betaInv(max(size(alpha), size(beta)))
        real(RKC)   :: betaInvRef(size(betaInv))
        real(RKC)   :: betaInc(NP)
        integer     :: fileUnit, i, j

        call setLinSpace(betaInc, 0._RKC, 1._RKC, fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "getBetaInv.RK.txt")
        do i = 1, NP
            betaInv = getBetaInv(betaInc(i), alpha, beta)
            write(fileUnit, "(*(g0,:,','))" ) betaInc(i), betaInv
        end do
        close(fileUnit)

        block
            use pm_kind, only: RKC => RK32
            character(*), parameter :: RKSTR = "RK32"
            real(RKC) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "getBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                betaInv = getBetaInv(real(betaInc(i), RKC), real(alpha, RKC), real(beta, RKC))
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKC => RK64
            character(*), parameter :: RKSTR = "RK64"
            real(RKC) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "getBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                betaInv = getBetaInv(real(betaInc(i), RKC), real(alpha, RKC), real(beta, RKC))
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKC => RK128
            character(*), parameter :: RKSTR = "RK128"
            real(RKC) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "getBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                betaInv = getBetaInv(real(betaInc(i), RKC), real(alpha, RKC), real(beta, RKC))
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

    end block

end program example