program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distPois, only: setPoisLogPMF

    implicit none

    real(RKC) :: logPMF(3)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setPoisLogPMF(logPMF(1), 0_IK, lambda = 2._RKC)")
                    call setPoisLogPMF(logPMF(1), 0_IK, lambda = 2._RKC)
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisLogPMF(logPMF(1), 0_IK, lambda = 2._RKC, logLambda = log(2._RKC))")
                    call setPoisLogPMF(logPMF(1), 0_IK, lambda = 2._RKC, logLambda = log(2._RKC))
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisLogPMF(logPMF(1), 1_IK, lambda = 2._RKC)")
                    call setPoisLogPMF(logPMF(1), 1_IK, lambda = 2._RKC)
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisLogPMF(logPMF(1), 2_IK, lambda = 2._RKC)")
                    call setPoisLogPMF(logPMF(1), 2_IK, lambda = 2._RKC)
    call disp%show("logPMF(1)")
    call disp%show( logPMF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisLogPMF(logPMF(1:3), [integer(IK) :: 0, 1, 2], lambda = 2._RKC)")
                    call setPoisLogPMF(logPMF(1:3), [integer(IK) :: 0, 1, 2], lambda = 2._RKC)
    call disp%show("logPMF(1:3)")
    call disp%show( logPMF(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setPoisLogPMF(logPMF(1:3), [integer(IK) :: 0, 1, 2], lambda = [0.1_RKC, 1._RKC, 10._RKC])")
                    call setPoisLogPMF(logPMF(1:3), [integer(IK) :: 0, 1, 2], lambda = [0.1_RKC, 1._RKC, 10._RKC])
    call disp%show("logPMF(1:3)")
    call disp%show( logPMF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        integer(IK), allocatable :: count(:)
        real(RKC) :: logPMF(4)
        integer :: fileUnit, i

        count = getRange(0_IK, 20_IK)
        open(newunit = fileUnit, file = "setPoisLogPMF.IK.txt")
        do i = 1, size(count)
            call setPoisLogPMF(logPMF, count(i), [.1_RKC, 1._RKC, 4._RKC, 10._RKC])
            write(fileUnit, "(*(g0,:,' '))" ) count(i), exp(logPMF)
        end do
        close(fileUnit)

    end block

end program example