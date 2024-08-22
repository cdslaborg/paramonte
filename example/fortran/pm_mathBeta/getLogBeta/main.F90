program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_mathBeta, only: getLogBeta

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized (lower) Incomplete Beta Function using its series representation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getLogBeta(alpha = 2._RK, beta = 3._RK) ! -2.484906649788000")
    call disp%show( getLogBeta(alpha = 2._RK, beta = 3._RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getLogBeta(alpha = 5._RK, beta = 3._RK) ! -4.653960350157524")
    call disp%show( getLogBeta(alpha = 5._RK, beta = 3._RK) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the regularized Incomplete Beta Function for a vector of input parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getLogBeta(alpha = [0.1_RK, 1._RK, 10._RK], beta = 3._RK) ! 2.158484749020289  -1.098612288668110  -6.492239835020470")
    call disp%show( getLogBeta(alpha = [0.1_RK, 1._RK, 10._RK], beta = 3._RK) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Incomplete Beta function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: getLinSpace
        integer(IK) , parameter :: NP = 1000_IK
        integer                 :: fileUnit, i
        real(RK)                :: A_RK(NP)

        A_RK = exp(getLinSpace(-2._RK, 2._RK, NP))
        open(newunit = fileUnit, file = "getLogBeta.RK.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))") A_RK(i), exp(getLogBeta(A_RK(i), [0.1_RK, 1.0_RK, 10._RK]))
        end do
        close(fileUnit)

    end block

end program example