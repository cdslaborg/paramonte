program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK128 ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_batse, only: getLogPbol

    implicit none

    integer(IK) :: info
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block

        use pm_kind, only: RKC => RK64

        call disp%skip()
        call disp%show("getLogPbol(logEpk = 2._RKC, logPF53 = 0._RKC)")
        call disp%show( getLogPbol(logEpk = 2._RKC, logPF53 = 0._RKC) )
        call disp%skip()

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        use pm_distBand, only: getBandEbreak, setBandEnergy
        real(RKC) :: logEpk(1000), energy, alpha = -1.1, beta = -2.3
        integer :: fileUnit, i, info

        call setLinSpace(logEpk, log(1.e-3_RKC), log(1.e5_RKC))
        open(newunit = fileUnit, file = "getLogPbol.RK.txt")
        do i = 1, size(logEpk)
            call setBandEnergy(energy, 0.001_RKC, 20000._RKC, 1._RKC, 50._RKC, 300._RKC, alpha, beta, getBandEbreak(alpha, beta, exp(logEpk(i))), info)
            if (info < 0) error stop
            write(fileUnit, "(*(g0,:,' '))" ) exp(logEpk(i)), exp(getLogPbol(logEpk(i), 0._RKC)), energy / 1.e9
        end do
        close(fileUnit)

    end block

end program example