program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distGeomCyclic, only: setGeomCyclicLogPMF

    implicit none

    real(RKC) :: logPMF(3)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setGeomCyclicLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKC), period = 1_IK)")
                    call setGeomCyclicLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKC), period = 1_IK)
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKC), period = 2_IK)")
                    call setGeomCyclicLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKC), period = 2_IK)
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKC), period = 20_IK)")
                    call setGeomCyclicLogPMF(logPMF(1), 1_IK, logProbSuccess = log(.2_RKC), period = 20_IK)
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogPMF(logPMF(1), 2_IK, logProbSuccess = log(.2_RKC), period = 2_IK)")
                    call setGeomCyclicLogPMF(logPMF(1), 2_IK, logProbSuccess = log(.2_RKC), period = 2_IK)
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogPMF(logPMF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log(.2_RKC), logProbFailure = log(1._RKC - .2_RKC), period = 4_IK)")
                    call setGeomCyclicLogPMF(logPMF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log(.2_RKC), logProbFailure = log(1._RKC - .2_RKC), period = 4_IK)
    call disp%show("logPMF(1:3)")
    call disp%show( logPMF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setGeomCyclicLogPMF(logPMF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log([0.1_RKC, .2_RKC, 1._RKC]), period = 5_IK)")
                    call setGeomCyclicLogPMF(logPMF(1:3), [integer(IK) :: 1, 2, 3], logProbSuccess = log([0.1_RKC, .2_RKC, 1._RKC]), period = 5_IK)
    call disp%show("logPMF(1:3)")
    call disp%show( logPMF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        integer(IK), allocatable :: stepSuccess(:)
        real(RKC) :: logPMF(4)
        integer :: fileUnit, i

        stepSuccess = getRange(1_IK, 10_IK)
        open(newunit = fileUnit, file = "setGeomCyclicLogPMF.IK.txt")
        do i = 1, size(stepSuccess)
            call setGeomCyclicLogPMF(logPMF, stepSuccess(i), log([.05_RKC, .25_RKC, .05_RKC, .25_RKC]), period = [10_IK, 10_IK, 10000_IK, 10000_IK])
            write(fileUnit, "(*(g0,:,','))" ) stepSuccess(i), exp(logPMF)
        end do
        close(fileUnit)

    end block

end program example