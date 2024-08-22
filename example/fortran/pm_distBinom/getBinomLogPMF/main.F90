program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: TKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBinom, only: getBinomLogPMF

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getBinomLogPMF(nsuc = 0_IK, ntry = [integer(IK) :: 0, 1, 2, 3, 10, 30], psuc = 0.5_TKG)")
    call disp%show( getBinomLogPMF(nsuc = 0_IK, ntry = [integer(IK) :: 0, 1, 2, 3, 10, 30], psuc = 0.5_TKG) )
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
        open(newunit = fileUnit, file = "getBinomLogPMF.IK.txt")
        do i = 1, size(nsuc)
            logPMF = getBinomLogPMF(nsuc(i), ntry, psuc)
            write(fileUnit, "(*(g0,:,' '))") nsuc(i), exp(logPMF)
        end do
        close(fileUnit)
    end block

end program example