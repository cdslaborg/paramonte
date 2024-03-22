program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK ! all processor kinds are supported.
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
    real(RKC)   :: alpha, beta, betaInc, betaInv
    real(RKC), allocatable :: alphas(:), betas(:), betaIncs(:), betaInvs(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    signed = .false.

    call disp%skip()
    call disp%show("signed")
    call disp%show( signed )
    call disp%show("betaInc = 0.5_RKC")
                    betaInc = 0.5_RKC
    call disp%show("alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)")
                    alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)
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
    call disp%show("betaInc = 0.5_RKC")
                    betaInc = 0.5_RKC
    call disp%show("alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)")
                    alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)
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
    call disp%show("betaInc = 1._RKC - epsilon(1._RKC)")
                    betaInc = 1._RKC - epsilon(1._RKC)
    call disp%show("alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)")
                    alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)
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
    call disp%show("betaIncs = [0._RKC, .5_RKC, 1._RKC]")
                    betaIncs = [0._RKC, .5_RKC, 1._RKC]
    call disp%show("alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)")
                    alpha = getUnifRand(0.1_RKC, 10._RKC); beta = getUnifRand(0.1_RKC, 10._RKC)
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
    call disp%show("betaIncs = [0._RKC, .5_RKC, 1._RKC]")
                    betaIncs = [0._RKC, .5_RKC, 1._RKC]
    call disp%show("alphas = [.1_RKC, 1._RKC, 10._RKC]")
                    alphas = [.1_RKC, 1._RKC, 10._RKC]
    call disp%show("beta = 3._RKC")
                    beta = 3._RKC
    call disp%show("call setResized(betaInvs, size(betaIncs, 1, IK))")
                    call setResized(betaInvs, size(betaIncs, 1, IK))
    call disp%show("call setBetaInv(betaInvs, betaIncs, alphas, beta, getLogBeta(alpha, beta), signed, infos)")
                    call setBetaInv(betaInvs, betaIncs, alphas, beta, getLogBeta(alpha, beta), signed, infos)
    call disp%show("if (any(infos /= 0)) error stop 'Beta inversion failed.'")
                    if (any(infos /= 0)) error stop 'Beta inversion failed.'
    call disp%show("betaInvs")
    call disp%show( betaInvs )
    call disp%show("getBetaInc(betaInvs, alpha = [.1_RKC, 1._RKC, 10._RKC], beta = 3._RKC, signed = signed)")
    call disp%show( getBetaInc(betaInvs, alpha = [.1_RKC, 1._RKC, 10._RKC], beta = 3._RKC, signed = signed) )
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
        integer(IK) :: infos(size(betaInv))
        integer     :: fileUnit, i, j

        call setLinSpace(betaInc, 0._RKC, 1._RKC, fopen = .true._LK, lopen = .true._LK)
        open(newunit = fileUnit, file = "setBetaInv.RK.txt")
        do i = 1, NP
            call setBetaInv(betaInv, betaInc(i), alpha, beta, getLogBeta(alpha, beta), signed, infos)
            if (any(infos /= 0)) error stop "Beta inversion failed."
            write(fileUnit, "(*(g0,:,','))" ) betaInc(i), betaInv
        end do
        close(fileUnit)

        block
            use pm_kind, only: RKC => RKS
            character(*), parameter :: RKSTR = "RKS"
            real(RKC) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "setBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                call setBetaInv(betaInv, real(betaInc(i), RKC), real(alpha, RKC), real(beta, RKC), getLogBeta(real(alpha, RKC), real(beta, RKC)), signed, infos)
                if (any(infos /= 0)) error stop "Beta inversion failed."
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKC => RKD
            character(*), parameter :: RKSTR = "RKD"
            real(RKC) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "setBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                call setBetaInv(betaInv, real(betaInc(i), RKC), real(alpha, RKC), real(beta, RKC), getLogBeta(real(alpha, RKC), real(beta, RKC)), signed, infos)
                if (any(infos /= 0)) error stop "Beta inversion failed."
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKC => RKH
            character(*), parameter :: RKSTR = "RKH"
            real(RKC) :: betaInv(max(size(alpha), size(beta)))
            open(newunit = fileUnit, file = "setBetaInv."//RKSTR//".abserr.txt")
            do i = 1, NP
                call setBetaInv(betaInv, real(betaInc(i), RKC), real(alpha, RKC), real(beta, RKC), getLogBeta(real(alpha, RKC), real(beta, RKC)), signed, infos)
                if (any(infos /= 0)) error stop "Beta inversion failed."
                betaInvRef = getBetaInv(betaInc(i), alpha, beta, signed = .true._LK)
                write(fileUnit, "(*(g0,:,','))") betaInc(i), abs(betaInv - merge(1 + betaInvRef, betaInvRef, betaInvRef < 0))
            end do
            close(fileUnit)
        end block

    end block

end program example