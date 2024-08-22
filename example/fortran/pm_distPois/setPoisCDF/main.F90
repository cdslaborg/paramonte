program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distPois, only: setPoisCDF
    use pm_err, only: setAsserted

    implicit none

    real(RKG) :: cdf(3)
    integer(IK) :: info(3)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1), 1._RKG, 2._RKG, info(1))")
                    call setPoisCDF(cdf(1), 1._RKG, 2._RKG, info(1))
    call disp%show("call setAsserted(.not. info(1) < 0)")
                    call setAsserted(.not. info(1) < 0)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1), 2._RKG, 2._RKG, info(1))")
                    call setPoisCDF(cdf(1), 2._RKG, 2._RKG, info(1))
    call disp%show("call setAsserted(.not. info(1) < 0)")
                    call setAsserted(.not. info(1) < 0)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1), 3._RKG, 2._RKG, info(1))")
                    call setPoisCDF(cdf(1), 3._RKG, 2._RKG, info(1))
    call disp%show("call setAsserted(.not. info(1) < 0)")
                    call setAsserted(.not. info(1) < 0)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1:3), 1._RKG + [real(RKG) :: 0, 1, 2], 2._RKG, info)")
                    call setPoisCDF(cdf(1:3), 1._RKG + [real(RKG) :: 0, 1, 2], 2._RKG, info)
    call disp%show("call setAsserted(.not. any(info < 0))")
                    call setAsserted(.not. any(info < 0))
    call disp%show("cdf(1:3)")
    call disp%show( cdf(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1:3), [real(RKG) :: 1, 2, 3], [0.1_RKG, 1._RKG, 10._RKG], info)")
                    call setPoisCDF(cdf(1:3), [real(RKG) :: 1, 2, 3], [0.1_RKG, 1._RKG, 10._RKG], info)
    call disp%show("call setAsserted(.not. any(info < 0))")
                    call setAsserted(.not. any(info < 0))
    call disp%show("cdf(1:3)")
    call disp%show( cdf(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: getLinSpace
        real(RKG), allocatable :: countP1(:)
        integer(IK) :: info(4)
        integer :: fileUnit, i
        real(RKG) :: cdf(4)

        countP1 = getLinSpace(1._RKG, 21._RKG, 21_IK)
        open(newunit = fileUnit, file = "setPoisCDF.IK.txt")
        do i = 1, size(countP1)
            call setPoisCDF(cdf, countP1(i), [.1_RKG, 1._RKG, 4._RKG, 10._RKG], info)
            call setAsserted(.not. any(info < 0))
            write(fileUnit, "(*(g0,:,' '))") countP1(i) - 1._RKG, cdf
        end do
        close(fileUnit)

    end block

end program example