program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distGeom, only: setGeomLogPMF

    implicit none

    real(RKG) :: logPMF(3)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setGeomLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKG))")
                    call setGeomLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKG))
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKG))")
                    call setGeomLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKG))
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKG))")
                    call setGeomLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKG))
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomLogPMF(logPMF(1), 2_IK, logProbSuccess = log(.2_RKG))")
                    call setGeomLogPMF(logPMF(1), 2_IK, logProbSuccess = log(.2_RKG))
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomLogPMF(logPMF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log(.2_RKG), logProbFailure = log(1._RKG - .2_RKG))")
                    call setGeomLogPMF(logPMF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log(.2_RKG), logProbFailure = log(1._RKG - .2_RKG))
    call disp%show("logPMF(1:3)")
    call disp%show( logPMF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomLogPMF(logPMF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log([0.1_RKG, .2_RKG, 1._RKG]))")
                    call setGeomLogPMF(logPMF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log([0.1_RKG, .2_RKG, 1._RKG]))
    call disp%show("logPMF(1:3)")
    call disp%show( logPMF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        integer(IK), allocatable :: stepSuccess(:)
        real(RKG) :: logPMF(4)
        integer :: fileUnit, i

        stepSuccess = getRange(1_IK, 10_IK)
        open(newunit = fileUnit, file = "setGeomLogPMF.IK.txt")
        do i = 1, size(stepSuccess)
            call setGeomLogPMF(logPMF, stepSuccess(i), log([.1_RKG, .2_RKG, .5_RKG, .8_RKG]))
            write(fileUnit, "(*(g0,:,' '))") stepSuccess(i), exp(logPMF)
        end do
        close(fileUnit)

    end block

end program example