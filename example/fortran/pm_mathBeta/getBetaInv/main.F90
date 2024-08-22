program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_mathBeta, only: getBetaInv
    use pm_distUnif, only: getUnifRand
    use pm_mathBeta, only: getBetaInc

    implicit none

    real(RKG), allocatable :: alpha, beta, betaInc(:), betaInv(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("betaInc = [0.5_RKG]")
                    betaInc = [0.5_RKG]
    call disp%show("alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)")
                    alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)
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
    call disp%show("betaInc = [0.5_RKG]")
                    betaInc = [0.5_RKG]
    call disp%show("alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)")
                    alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)
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
    call disp%show("betaInc = [1._RKG - epsilon(1._RKG)]")
                    betaInc = [1._RKG - epsilon(1._RKG)]
    call disp%show("alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)")
                    alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)
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
    call disp%show("betaInc = [0._RKG, .5_RKG, 1._RKG]")
                    betaInc = [0._RKG, .5_RKG, 1._RKG]
    call disp%show("alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)")
                    alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)
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
    call disp%show("betaInc = [0._RKG, .5_RKG, 1._RKG]")
                    betaInc = [0._RKG, .5_RKG, 1._RKG]
    call disp%show("betaInv = getBetaInv(betaInc, alpha = [.1_RKG, 1._RKG, 10._RKG], beta = 3._RKG)")
                    betaInv = getBetaInv(betaInc, alpha = [.1_RKG, 1._RKG, 10._RKG], beta = 3._RKG)
    call disp%show("betaInv")
    call disp%show( betaInv )
    call disp%show("getBetaInc(betaInv, alpha = [.1_RKG, 1._RKG, 10._RKG], beta = 3._RKG)")
    call disp%show( getBetaInc(betaInv, alpha = [.1_RKG, 1._RKG, 10._RKG], beta = 3._RKG) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_kind, only: RKG => RKH, RKH
        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 1000
        real(RKG)   , parameter :: alpha(*) = [0.1_RKG, 10._RKG, 1._RKG, 0.1_RKG, 10._RKG], beta(*) = [0.1_RKG, 0.1_RKG, 1._RKG, 10._RKG, 10._RKG]
        real(RKG)   :: betaInv(max(size(alpha), size(beta)))
        real(RKG)   :: betaInvRef(size(betaInv))
        real(RKG)   :: betaInc(NP)
        integer     :: fileUnit, i, j

        call setLinSpace(betaInc, 0._RKG, 1._RKG, fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "getBetaInv.RK.txt")
        do i = 1, NP
            betaInv = getBetaInv(betaInc(i), alpha, beta)
            write(fileUnit, "(*(g0,:,','))") betaInc(i), betaInv
        end do
        close(fileUnit)

        block
            use pm_kind, only: RKG => RKS
            character(*), parameter :: RKSTR = "RKS"
            real(RKG) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "getBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                betaInv = getBetaInv(real(betaInc(i), RKG), real(alpha, RKG), real(beta, RKG))
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKG => RKD
            character(*), parameter :: RKSTR = "RKD"
            real(RKG) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "getBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                betaInv = getBetaInv(real(betaInc(i), RKG), real(alpha, RKG), real(beta, RKG))
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKG => RKH
            character(*), parameter :: RKSTR = "RKH"
            real(RKG) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "getBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                betaInv = getBetaInv(real(betaInc(i), RKG), real(alpha, RKG), real(beta, RKG))
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

    end block

end program example