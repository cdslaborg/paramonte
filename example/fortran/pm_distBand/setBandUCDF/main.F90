program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBand, only: setBandUCDF

    implicit none

    real(RKG) :: ucdf(4)
    integer(IK) :: info(4)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setBandUCDF(ucdf(1), lb = .01_RKG, ub = 10._RKG, alpha = +2._RKG, beta = -3._RKG, ebreak = 1._RKG, info = info(1))")
                    call setBandUCDF(ucdf(1), lb = .01_RKG, ub = 10._RKG, alpha = +2._RKG, beta = -3._RKG, ebreak = 1._RKG, info = info(1))
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
        real(RKG) :: ucdf(4), ub(NP)
        integer :: fileUnit, i

        call setLinSpace(ub, 0.01_RKG, 10._RKG)
        open(newunit = fileUnit, file = "setBandUCDF.RK.txt")
        do i = 1, NP
            call setBandUCDF(ucdf, 0.01_RKG, ub(i), [.5_RKG, 1.5_RKG, -.5_RKG, -1.1_RKG], -[.5_RKG, 1._RKG, 2._RKG, 3._RKG], [.5_RKG, 1.0_RKG, 2.0_RKG, 5._RKG], info = info)
            if (any(info < 0)) error stop
            write(fileUnit, "(*(g0,:,' '))") ub(i), ucdf
        end do
        close(fileUnit)

    end block

end program example