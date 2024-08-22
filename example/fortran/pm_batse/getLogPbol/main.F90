program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_batse, only: getLogPbol

    implicit none

    integer(IK) :: info
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block

        use pm_kind, only: RKG => RKD

        call disp%skip()
        call disp%show("getLogPbol(logEpk = 2._RKG, logPF53 = 0._RKG)")
        call disp%show( getLogPbol(logEpk = 2._RKG, logPF53 = 0._RKG) )
        call disp%skip()

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        use pm_distBand, only: getBandEbreak, setBandEnergy
        real(RKG) :: logEpk(1000), energy, alpha = -1.1, beta = -2.3
        integer :: fileUnit, i, info

        call setLinSpace(logEpk, log(1.e-3_RKG), log(1.e5_RKG))
        open(newunit = fileUnit, file = "getLogPbol.RK.txt")
        do i = 1, size(logEpk)
            call setBandEnergy(energy, 0.001_RKG, 20000._RKG, 1._RKG, 50._RKG, 300._RKG, alpha, beta, getBandEbreak(alpha, beta, exp(logEpk(i))), info)
            if (info < 0) error stop
            write(fileUnit, "(*(g0,:,' '))") exp(logEpk(i)), exp(getLogPbol(logEpk(i), 0._RKG)), energy / 1.e9
        end do
        close(fileUnit)

    end block

end program example