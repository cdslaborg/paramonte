program example

    use pm_kind, only: SK, IK
    use pm_arraySpace, only: getLinSpace
    use pm_distGeomCyclic, only: getGeomCyclicRand
    use pm_distGeomCyclic, only: isFailedGeomCyclicFit
    use pm_distGeomCyclic, only: getGeomCyclicLogPMF
    use pm_arrayUnique, only: setUnique
    use pm_io, only: getErrTableWrite
    use pm_io, only: display_type

    implicit none

    integer(IK) :: period, i
    type(display_type) :: disp
    integer(IK), allocatable :: rand(:)
    integer(IK), allocatable :: freqSuccess(:)
    integer(IK), allocatable :: stepSuccess(:)
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: RKG => RKH
        real(RKG) :: probSuccess, probSuccessFit, normFacFit
        real(RKG), allocatable :: x(:), y(:)
        call disp%skip()
        call disp%show("probSuccess = .1; period = 20")
                        probSuccess = .1; period = 20
        call disp%show("rand = getGeomCyclicRand(probSuccess, [(period, i = 1, 1000)])")
                        rand = getGeomCyclicRand(probSuccess, [(period, i = 1, 1000)])
        call disp%show("if (0 /= getErrTableWrite(file = 'isFailedGeomCyclicFit.rnd.txt', table = rand)) error stop 'table output failed.'")
                        if (0 /= getErrTableWrite(file = 'isFailedGeomCyclicFit.rnd.txt', table = rand)) error stop 'table output failed.'
        call disp%show("call setUnique(getGeomCyclicRand(probSuccess, [(period, i = 1, 1000)]), unique = stepSuccess, count = freqSuccess, order = -1_IK)")
                        call setUnique(getGeomCyclicRand(probSuccess, [(period, i = 1, 1000)]), unique = stepSuccess, count = freqSuccess, order = -1_IK)
        call disp%show("stepSuccess")
        call disp%show( stepSuccess , format = SK_"(sp,20(g0,:,', '))")
        call disp%show("freqSuccess")
        call disp%show( freqSuccess , format = SK_"(sp,20(g0,:,', '))")
        call disp%show("if (isFailedGeomCyclicFit(stepSuccess, freqSuccess, period, probSuccessFit, normFacFit)) error stop 'Fitting failed.'")
                        if (isFailedGeomCyclicFit(stepSuccess, freqSuccess, period, probSuccessFit, normFacFit)) error stop 'Fitting failed.'
        call disp%show("[probSuccessFit, normFacFit]")
        call disp%show( [probSuccessFit, normFacFit] )
        call disp%show("x = getLinSpace(1._RKG, real(period, RKG), 500_IK) ! for visualization.")
                        x = getLinSpace(1._RKG, real(period, RKG), 500_IK) ! for visualization.
        call disp%show("y = normFacFit * probSuccessFit * (1 - probSuccessFit)**(x - 1) / (1 - (1 - probSuccessFit)**period) ! for visualization.")
                        y = normFacFit * probSuccessFit * (1 - probSuccessFit)**(x - 1) / (1 - (1 - probSuccessFit)**period) ! for visualization.
        call disp%show("y = normFacFit * exp(getGeomCyclicLogPMF(x, probSuccessFit, period)) ! for visualization.")
                        y = normFacFit * exp(getGeomCyclicLogPMF(x, probSuccessFit, period)) ! for visualization.
        call disp%show("if (0 /= getErrTableWrite(file = 'isFailedGeomCyclicFit.fit.txt', table = reshape([x, y], [size(x), 2])))    error stop 'table output failed.'")
                        if (0 /= getErrTableWrite(file = 'isFailedGeomCyclicFit.fit.txt', table = reshape([x, y], [size(x), 2])))    error stop 'table output failed.'
        call disp%skip()
    end block

end program example