program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK128 ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBand, only: setBandUCDF

    implicit none

    real(RKC) :: ucdf(4)
    integer(IK) :: info(4)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setBandUCDF(ucdf(1), lb = .01_RKC, ub = 10._RKC, alpha = +2._RKC, beta = -3._RKC, ebreak = 1._RKC, info = info(1))")
                    call setBandUCDF(ucdf(1), lb = .01_RKC, ub = 10._RKC, alpha = +2._RKC, beta = -3._RKC, ebreak = 1._RKC, info = info(1))
    call disp%show("if (info(1) < 0) error stop")
                    if (info(1) < 0) error stop
    call disp%show("ucdf(1)")
    call disp%show( ucdf(1) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 1000_IK
        real(RKC) :: ucdf(4), ub(NP)
        integer :: fileUnit, i

        call setLinSpace(ub, 0.01_RKC, 10._RKC)
        open(newunit = fileUnit, file = "setBandUCDF.RK.txt")
        do i = 1, NP
            call setBandUCDF(ucdf, 0.01_RKC, ub(i), [.5_RKC, 1.5_RKC, -.5_RKC, -1.1_RKC], -[.5_RKC, 1._RKC, 2._RKC, 3._RKC], [.5_RKC, 1.0_RKC, 2.0_RKC, 5._RKC], info = info)
            if (any(info < 0)) error stop
            write(fileUnit, "(*(g0,:,' '))" ) ub(i), ucdf
        end do
        close(fileUnit)

    end block

end program example