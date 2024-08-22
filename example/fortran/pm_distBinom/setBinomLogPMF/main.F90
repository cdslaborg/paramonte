program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: TKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBinom, only: setBinomLogPMF

    implicit none

    type(display_type) :: disp
    real(TKG), allocatable :: logPMF(:)
    integer(IK), parameter :: ntry(*) = [integer(IK) :: 0, 1, 2, 3, 10, 30]
    disp = display_type(file = "main.out.F90")


    call disp%skip()
    call disp%show("ntry")
    call disp%show( ntry )
    call disp%show("allocate(logPMF(size(ntry)))")
                    allocate(logPMF(size(ntry)))
    call disp%show("call setBinomLogPMF(logPMF, nsuc = 0_IK, ntry = ntry, logp = log(0.5_TKG), logq = log(0.5_TKG))")
                    call setBinomLogPMF(logPMF, nsuc = 0_IK, ntry = ntry, logp = log(0.5_TKG), logq = log(0.5_TKG))
    call disp%show("logPMF")
    call disp%show( logPMF )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_except, only: getNAN
        use pm_arrayRange, only: getRange
        integer(IK), parameter :: ntry = 60
        real(TKG), parameter :: psuc(*) = [real(TKG) :: .1, .5, .9]
        integer(IK), allocatable :: nsuc(:)
        real(TKG) :: logPMF(size(psuc))
        integer(IK) :: fileUnit, i
        nsuc = getRange(0_IK, ntry)
        open(newunit = fileUnit, file = "setBinomLogPMF.IK.txt")
        do i = 1, size(nsuc)
            call setBinomLogPMF(logPMF, nsuc(i), ntry, log(psuc), log(1 - psuc))
            write(fileUnit, "(*(g0,:,' '))") nsuc(i), exp(logPMF)
        end do
        close(fileUnit)
    end block

end program example