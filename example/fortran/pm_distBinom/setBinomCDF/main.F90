program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: TKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBinom, only: setBinomCDF

    implicit none

    type(display_type) :: disp
    real(TKG), allocatable :: cdf(:)
    integer(IK), allocatable :: info(:)
    integer(IK), parameter :: ntry(*) = [integer(IK) :: 0, 1, 2, 3, 10, 30]
    disp = display_type(file = "main.out.F90")


    call disp%skip()
    call disp%show("ntry")
    call disp%show( ntry )
    call disp%show("allocate(cdf(size(ntry)), info(size(ntry)))")
                    allocate(cdf(size(ntry)), info(size(ntry)))
    call disp%show("call setBinomCDF(cdf, nsuc = 0_IK, ntry = ntry, psuc = 0.5_TKG, info = info)")
                    call setBinomCDF(cdf, nsuc = 0_IK, ntry = ntry, psuc = 0.5_TKG, info = info)
    call disp%show("if (any(info /= 0)) error stop 'setBinomCDF() failed.'")
                    if (any(info /= 0)) error stop 'setBinomCDF() failed.' 
    call disp%show("info")
    call disp%show( info )
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arrayRange, only: getRange
        integer(IK), parameter :: ntry = 60
        real(TKG), parameter :: psuc(*) = [real(TKG) :: .1, .5, .9]
        integer(IK), allocatable :: nsuc(:), info(:)
        real(TKG) :: cdf(size(psuc))
        integer(IK) :: fileUnit, i
        nsuc = getRange(0_IK, ntry)
        allocate(info(size(psuc)))
        open(newunit = fileUnit, file = "setBinomCDF.IK.txt")
        do i = 1, size(nsuc)
            call setBinomCDF(cdf, nsuc(i), ntry, psuc, info)
            if (any(info /= 0)) error stop "setBinomCDF() Failed."
            write(fileUnit, "(*(g0,:,' '))") nsuc(i), cdf
        end do
        close(fileUnit)
    end block

end program example