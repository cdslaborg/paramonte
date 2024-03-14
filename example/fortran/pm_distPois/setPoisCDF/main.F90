program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distPois, only: setPoisCDF
    use pm_err, only: setAsserted

    implicit none

    real(RKC) :: cdf(3)
    integer(IK) :: info(3)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1), 1._RKC, log_gamma(1._RKC), 2._RKC, info(1))")
                    call setPoisCDF(cdf(1), 1._RKC, log_gamma(1._RKC), 2._RKC, info(1))
    call disp%show("call setAsserted(.not. info(1) < 0)")
                    call setAsserted(.not. info(1) < 0)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1), 2._RKC, log_gamma(2._RKC), 2._RKC, info(1))")
                    call setPoisCDF(cdf(1), 2._RKC, log_gamma(2._RKC), 2._RKC, info(1))
    call disp%show("call setAsserted(.not. info(1) < 0)")
                    call setAsserted(.not. info(1) < 0)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1), 3._RKC, log_gamma(3._RKC), 2._RKC, info(1))")
                    call setPoisCDF(cdf(1), 3._RKC, log_gamma(3._RKC), 2._RKC, info(1))
    call disp%show("call setAsserted(.not. info(1) < 0)")
                    call setAsserted(.not. info(1) < 0)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1:3), 1._RKC + [real(RKC) :: 0, 1, 2], log_gamma(1._RKC + [real(RKC) :: 0, 1, 2]), 2._RKC, info)")
                    call setPoisCDF(cdf(1:3), 1._RKC + [real(RKC) :: 0, 1, 2], log_gamma(1._RKC + [real(RKC) :: 0, 1, 2]), 2._RKC, info)
    call disp%show("call setAsserted(.not. any(info < 0))")
                    call setAsserted(.not. any(info < 0))
    call disp%show("cdf(1:3)")
    call disp%show( cdf(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisCDF(cdf(1:3), [real(RKC) :: 1, 2, 3], log_gamma([real(RKC) :: 1, 2, 3]), [0.1_RKC, 1._RKC, 10._RKC], info)")
                    call setPoisCDF(cdf(1:3), [real(RKC) :: 1, 2, 3], log_gamma([real(RKC) :: 1, 2, 3]), [0.1_RKC, 1._RKC, 10._RKC], info)
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
        real(RKC), allocatable :: countP1(:)
        integer(IK) :: info(4)
        integer :: fileUnit, i
        real(RKC) :: cdf(4)

        countP1 = getLinSpace(1._RKC, 21._RKC, 21_IK)
        open(newunit = fileUnit, file = "setPoisCDF.IK.txt")
        do i = 1, size(countP1)
            call setPoisCDF(cdf, countP1(i), log_gamma(countP1(i)), [.1_RKC, 1._RKC, 4._RKC, 10._RKC], info)
            call setAsserted(.not. any(info < 0))
            write(fileUnit, "(*(g0,:,' '))" ) countP1(i) - 1._RKC, cdf
        end do
        close(fileUnit)

    end block

end program example