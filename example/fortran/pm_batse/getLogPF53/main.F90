program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_batse, only: getLogPF53

    implicit none

    integer(IK) :: info
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block

        use pm_kind, only: RKG => RKD

        call disp%skip()
        call disp%show("getLogPF53(logEpk = 2._RKG, logPbol = 0._RKG)")
        call disp%show( getLogPF53(logEpk = 2._RKG, logPbol = 0._RKG) )
        call disp%skip()

    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        use pm_distBand, only: getBandEbreak, setBandPhoton
        real(RKG) :: logEpk(1000), photon, alpha = -1.1, beta = -2.3
        integer :: fileUnit, i

        call setLinSpace(logEpk, log(1.e-3_RKG), log(1.e5_RKG))
        open(newunit = fileUnit, file = "getLogPF53.RK.txt")
        do i = 1, size(logEpk)
            call setBandPhoton(photon, 50._RKG, 300._RKG, 1._RKG, 0.001_RKG, 20000._RKG, alpha, beta, getBandEbreak(alpha, beta, exp(logEpk(i))), info)
            if (info < 0) error stop
            write(fileUnit, "(*(g0,:,' '))") exp(logEpk(i)), exp(getLogPF53(logEpk(i), log(1.e-6_RKG))), photon
        end do
        close(fileUnit)

    end block

end program example