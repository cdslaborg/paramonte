program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_mathBeta, only: getBetaInv
    use pm_mathBeta, only: setBetaInv
    use pm_distUnif, only: getUnifRand
    use pm_mathBeta, only: getBetaInc
    use pm_mathBeta, only: getLogBeta
    use pm_arrayResize, only: setResized

    implicit none

    logical(LK) :: signed
    integer(IK) :: info, infos(3)
    real(RKG)   :: alpha, beta, betaInc, betaInv
    real(RKG), allocatable :: alphas(:), betas(:), betaIncs(:), betaInvs(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    signed = .false.

    call disp%skip()
    call disp%show("signed")
    call disp%show( signed )
    call disp%show("betaInc = 0.5_RKG")
                    betaInc = 0.5_RKG
    call disp%show("alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)")
                    alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)
    call disp%show("[alpha, beta]")
    call disp%show( [alpha, beta] )
    call disp%show("call setBetaInv(betaInv, betaInc, alpha, beta, getLogBeta(alpha, beta), signed, info)")
                    call setBetaInv(betaInv, betaInc, alpha, beta, getLogBeta(alpha, beta), signed, info)
    call disp%show("if (info /= 0) error stop 'Beta inversion failed.'")
                    if (info /= 0) error stop 'Beta inversion failed.'
    call disp%show("betaInv")
    call disp%show( betaInv )
    call disp%show("getBetaInc(betaInv, alpha, beta, signed)")
    call disp%show( getBetaInc(betaInv, alpha, beta, signed) )
    call disp%skip()

    call disp%skip()
    call disp%show("signed")
    call disp%show( signed )
    call disp%show("betaInc = 0.5_RKG")
                    betaInc = 0.5_RKG
    call disp%show("alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)")
                    alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)
    call disp%show("[alpha, beta]")
    call disp%show( [alpha, beta] )
    call disp%show("call setBetaInv(betaInv, betaInc, alpha, beta, getLogBeta(alpha, beta), signed, info)")
                    call setBetaInv(betaInv, betaInc, alpha, beta, getLogBeta(alpha, beta), signed, info)
    call disp%show("if (info /= 0) error stop 'Beta inversion failed.'")
                    if (info /= 0) error stop 'Beta inversion failed.'
    call disp%show("betaInv")
    call disp%show( betaInv )
    call disp%show("getBetaInc(betaInv, alpha, beta, signed)")
    call disp%show( getBetaInc(betaInv, alpha, beta, signed) )
    call disp%skip()

    call disp%skip()
    call disp%show("signed")
    call disp%show( signed )
    call disp%show("betaInc = 1._RKG - epsilon(1._RKG)")
                    betaInc = 1._RKG - epsilon(1._RKG)
    call disp%show("alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)")
                    alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)
    call disp%show("[alpha, beta]")
    call disp%show( [alpha, beta] )
    call disp%show("call setBetaInv(betaInv, betaInc, alpha, beta, getLogBeta(alpha, beta), signed, info)")
                    call setBetaInv(betaInv, betaInc, alpha, beta, getLogBeta(alpha, beta), signed, info)
    call disp%show("if (info /= 0) error stop 'Beta inversion failed.'")
                    if (info /= 0) error stop 'Beta inversion failed.'
    call disp%show("betaInv")
    call disp%show( betaInv )
    call disp%show("getBetaInc(betaInv, alpha, beta, signed)")
    call disp%show( getBetaInc(betaInv, alpha, beta, signed) )
    call disp%skip()

    call disp%skip()
    call disp%show("signed")
    call disp%show( signed )
    call disp%show("betaIncs = [0._RKG, .5_RKG, 1._RKG]")
                    betaIncs = [0._RKG, .5_RKG, 1._RKG]
    call disp%show("alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)")
                    alpha = getUnifRand(0.1_RKG, 10._RKG); beta = getUnifRand(0.1_RKG, 10._RKG)
    call disp%show("[alpha, beta]")
    call disp%show( [alpha, beta] )
    call disp%show("call setResized(betaInvs, size(betaIncs, 1, IK))")
                    call setResized(betaInvs, size(betaIncs, 1, IK))
    call disp%show("call setBetaInv(betaInvs, betaIncs, alpha, beta, getLogBeta(alpha, beta), signed, infos)")
                    call setBetaInv(betaInvs, betaIncs, alpha, beta, getLogBeta(alpha, beta), signed, infos)
    call disp%show("if (any(infos /= 0)) error stop 'Beta inversion failed.'")
                    if (any(infos /= 0)) error stop 'Beta inversion failed.'
    call disp%show("betaInvs")
    call disp%show( betaInvs )
    call disp%show("getBetaInc(betaInvs, alpha, beta, signed)")
    call disp%show( getBetaInc(betaInvs, alpha, beta, signed) )
    call disp%skip()

    call disp%skip()
    call disp%show("signed")
    call disp%show( signed )
    call disp%show("betaIncs = [0._RKG, .5_RKG, 1._RKG]")
                    betaIncs = [0._RKG, .5_RKG, 1._RKG]
    call disp%show("alphas = [.1_RKG, 1._RKG, 10._RKG]")
                    alphas = [.1_RKG, 1._RKG, 10._RKG]
    call disp%show("beta = 3._RKG")
                    beta = 3._RKG
    call disp%show("call setResized(betaInvs, size(betaIncs, 1, IK))")
                    call setResized(betaInvs, size(betaIncs, 1, IK))
    call disp%show("call setBetaInv(betaInvs, betaIncs, alphas, beta, getLogBeta(alpha, beta), signed, infos)")
                    call setBetaInv(betaInvs, betaIncs, alphas, beta, getLogBeta(alpha, beta), signed, infos)
    call disp%show("if (any(infos /= 0)) error stop 'Beta inversion failed.'")
                    if (any(infos /= 0)) error stop 'Beta inversion failed.'
    call disp%show("betaInvs")
    call disp%show( betaInvs )
    call disp%show("getBetaInc(betaInvs, alpha = [.1_RKG, 1._RKG, 10._RKG], beta = 3._RKG, signed = signed)")
    call disp%show( getBetaInc(betaInvs, alpha = [.1_RKG, 1._RKG, 10._RKG], beta = 3._RKG, signed = signed) )
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
        integer(IK) :: infos(size(betaInv))
        integer     :: fileUnit, i, j

        call setLinSpace(betaInc, 0._RKG, 1._RKG, fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "setBetaInv.RK.txt")
        do i = 1, NP
            call setBetaInv(betaInv, betaInc(i), alpha, beta, getLogBeta(alpha, beta), signed, infos)
            if (any(infos /= 0)) error stop "Beta inversion failed."
            write(fileUnit, "(*(g0,:,','))") betaInc(i), betaInv
        end do
        close(fileUnit)

        block
            use pm_kind, only: RKG => RKS
            character(*), parameter :: RKSTR = "RKS"
            real(RKG) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "setBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                call setBetaInv(betaInv, real(betaInc(i), RKG), real(alpha, RKG), real(beta, RKG), getLogBeta(real(alpha, RKG), real(beta, RKG)), signed, infos)
                if (any(infos /= 0)) error stop "Beta inversion failed."
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKG => RKD
            character(*), parameter :: RKSTR = "RKD"
            real(RKG) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "setBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                call setBetaInv(betaInv, real(betaInc(i), RKG), real(alpha, RKG), real(beta, RKG), getLogBeta(real(alpha, RKG), real(beta, RKG)), signed, infos)
                if (any(infos /= 0)) error stop "Beta inversion failed."
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKG => RKH
            character(*), parameter :: RKSTR = "RKH"
            real(RKG) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "setBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                call setBetaInv(betaInv, real(betaInc(i), RKG), real(alpha, RKG), real(beta, RKG), getLogBeta(real(alpha, RKG), real(beta, RKG)), signed, infos)
                if (any(infos /= 0)) error stop "Beta inversion failed."
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

    end block

end program example