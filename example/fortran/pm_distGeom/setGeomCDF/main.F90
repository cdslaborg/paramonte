program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distGeom, only: setGeomCDF

    implicit none

    real(RKG) :: cdf(3)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setGeomCDF(cdf(1), 1_IK, probSuccess = .2_RKG)")
                    call setGeomCDF(cdf(1), 1_IK, probSuccess = .2_RKG)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCDF(cdf(1), 1_IK, probSuccess = .2_RKG)")
                    call setGeomCDF(cdf(1), 1_IK, probSuccess = .2_RKG)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCDF(cdf(1), 1_IK, probSuccess = .2_RKG)")
                    call setGeomCDF(cdf(1), 1_IK, probSuccess = .2_RKG)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCDF(cdf(1), 2_IK, probSuccess = .2_RKG)")
                    call setGeomCDF(cdf(1), 2_IK, probSuccess = .2_RKG)
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCDF(cdf(1:3), [integer(IK) :: 1, 2, 3], probSuccess = .2_RKG)")
                    call setGeomCDF(cdf(1:3), [integer(IK) :: 1, 2, 3], probSuccess = .2_RKG)
    call disp%show("cdf(1:3)")
    call disp%show( cdf(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCDF(cdf(1:3), [integer(IK) :: 1, 2, 3], probSuccess = [0.1_RKG, .2_RKG, 1._RKG])")
                    call setGeomCDF(cdf(1:3), [integer(IK) :: 1, 2, 3], probSuccess = [0.1_RKG, .2_RKG, 1._RKG])
    call disp%show("cdf(1:3)")
    call disp%show( cdf(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        integer(IK), allocatable :: stepSuccess(:)
        real(RKG) :: cdf(4)
        integer :: fileUnit, i

        stepSuccess = getRange(1_IK, 10_IK)
        open(newunit = fileUnit, file = "setGeomCDF.IK.txt")
        do i = 1, size(stepSuccess)
            call setGeomCDF(cdf, stepSuccess(i), [.1_RKG, .2_RKG, .5_RKG, .8_RKG])
            write(fileUnit, "(*(g0,:,' '))") stepSuccess(i), cdf
        end do
        close(fileUnit)

    end block

end program example