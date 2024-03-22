program example

    use pm_kind, only: SK, IK, RKH ! All real kinds are supported.
    use pm_cosmology, only: getHubbleParamNormedSq
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getHubbleParamNormedSq(zplus1 = 1.1)")
    call disp%show( getHubbleParamNormedSq(zplus1 = 1.1) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHubbleParamNormedSq(zplus1 = 1.1d0)")
    call disp%show( getHubbleParamNormedSq(zplus1 = 1.1d0) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHubbleParamNormedSq(zplus1 = 1.1_RKH)")
    call disp%show( getHubbleParamNormedSq(zplus1 = 1.1_RKH) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHubbleParamNormedSq(zplus1 = [real :: 1, 2, 3, 4])")
    call disp%show( getHubbleParamNormedSq(zplus1 = [real :: 1, 2, 3, 4]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHubbleParamNormedSq(zplus1 = [real :: 1, 2, 3, 4], omegaM = 0.3, omegaL = 0.7)")
    call disp%show( getHubbleParamNormedSq(zplus1 = [real :: 1, 2, 3, 4], omegaM = 0.3, omegaL = 0.7) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHubbleParamNormedSq(zplus1 = [real :: 1, 2, 3, 4], omegaM = 0.3, omegaL = 0.5, omegaR = 0.2)")
    call disp%show( getHubbleParamNormedSq(zplus1 = [real :: 1, 2, 3, 4], omegaM = 0.3, omegaL = 0.5, omegaR = 0.2) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHubbleParamNormedSq(zplus1 = [real :: 1, 2, 3, 4], omegaM = 0.3, omegaL = 0.5, omegaR = 0.1, omegaK = 0.1)")
    call disp%show( getHubbleParamNormedSq(zplus1 = [real :: 1, 2, 3, 4], omegaM = 0.3, omegaL = 0.5, omegaR = 0.1, omegaK = 0.1) )
    call disp%skip()

    ! Visualize the Hubble Parameter.

    block

        use pm_cosmology, only:  HUBBLE_CONST
        use pm_arraySpace, only: getLogSpace

        integer(IK) , parameter :: NP = 1000_IK
        real                    :: zplus1(NP)
        integer(IK)             :: i, fileUnit

        open(newunit = fileUnit, file = "getHubbleParamNormedSq.RK.txt")
        zplus1 = getLogSpace(logx1 = log(1.), logx2 = log(10.), count = NP)
        do i = 1_IK, size(zplus1, 1, IK)
            write(fileUnit, "(*(g0,:,', '))") zplus1(i) &
                                            , real(HUBBLE_CONST) * sqrt(getHubbleParamNormedSq(zplus1(i))) &
                                            , real(HUBBLE_CONST) * sqrt(getHubbleParamNormedSq(zplus1(i), 1.0, 0.0)) &
                                            , real(HUBBLE_CONST) * sqrt(getHubbleParamNormedSq(zplus1(i), 0.0, 1.0)) &
                                            , real(HUBBLE_CONST) * sqrt(getHubbleParamNormedSq(zplus1(i), .30, .30, .40)) &
                                            , real(HUBBLE_CONST) * sqrt(getHubbleParamNormedSq(zplus1(i), .25, .25, .25, .25))
        end do
        close(fileUnit)

    end block

end program example