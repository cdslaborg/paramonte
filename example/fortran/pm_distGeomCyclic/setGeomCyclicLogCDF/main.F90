program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distGeomCyclic, only: setGeomCyclicLogCDF

    implicit none

    real(RKG) :: logCDF(3)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setGeomCyclicLogCDF(logCDF(1), 1_IK, logProbSuccess = log(.2_RKG), period = 1_IK)")
                    call setGeomCyclicLogCDF(logCDF(1), 1_IK, logProbSuccess = log(.2_RKG), period = 1_IK)
    call disp%show("logCDF(1)")
    call disp%show( logCDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogCDF(logCDF(1), 1_IK, logProbSuccess = log(.2_RKG), period = 2_IK)")
                    call setGeomCyclicLogCDF(logCDF(1), 1_IK, logProbSuccess = log(.2_RKG), period = 2_IK)
    call disp%show("logCDF(1)")
    call disp%show( logCDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogCDF(logCDF(1), 1_IK, logProbSuccess = log(.2_RKG), period = 20_IK)")
                    call setGeomCyclicLogCDF(logCDF(1), 1_IK, logProbSuccess = log(.2_RKG), period = 20_IK)
    call disp%show("logCDF(1)")
    call disp%show( logCDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogCDF(logCDF(1), 2_IK, logProbSuccess = log(.2_RKG), period = 2_IK)")
                    call setGeomCyclicLogCDF(logCDF(1), 2_IK, logProbSuccess = log(.2_RKG), period = 2_IK)
    call disp%show("logCDF(1)")
    call disp%show( logCDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogCDF(logCDF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log(.2_RKG), logProbFailure = log(1._RKG - .2_RKG), period = 4_IK)")
                    call setGeomCyclicLogCDF(logCDF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log(.2_RKG), logProbFailure = log(1._RKG - .2_RKG), period = 4_IK)
    call disp%show("logCDF(1:3)")
    call disp%show( logCDF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogCDF(logCDF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log([0.1_RKG, .2_RKG, 1._RKG]), period = 5_IK)")
                    call setGeomCyclicLogCDF(logCDF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log([0.1_RKG, .2_RKG, 1._RKG]), period = 5_IK)
    call disp%show("logCDF(1:3)")
    call disp%show( logCDF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        integer(IK), allocatable :: stepSuccess(:)
        real(RKG) :: logCDF(4)
        integer :: fileUnit, i

        stepSuccess = getRange(1_IK, 10_IK)
        open(newunit = fileUnit, file = "setGeomCyclicLogCDF.IK.txt")
        do i = 1, size(stepSuccess)
            call setGeomCyclicLogCDF(logCDF, stepSuccess(i), log([.05_RKG, .25_RKG, .05_RKG, .25_RKG]), period = [10_IK, 10_IK, 10000_IK, 10000_IK])
            write(fileUnit, "(*(g0,:,','))") stepSuccess(i), exp(logCDF)
        end do
        close(fileUnit)

    end block

end program example