program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBeta, only: getBetaCDF

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getBetaCDF(x = 0._RKC, alpha = 2._RKC, beta = 3._RKC)")
    call disp%show( getBetaCDF(x = 0._RKC, alpha = 2._RKC, beta = 3._RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaCDF(x = .5_RKC, alpha = 2._RKC, beta = 3._RKC)")
    call disp%show( getBetaCDF(x = .5_RKC, alpha = 2._RKC, beta = 3._RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaCDF(x = 1._RKC, alpha = 2._RKC, beta = 3._RKC)")
    call disp%show( getBetaCDF(x = 1._RKC, alpha = 2._RKC, beta = 3._RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaCDF(x = [0._RKC, 0.5_RKC, 1._RKC], alpha = 2._RKC, beta = 3._RKC)")
    call disp%show( getBetaCDF(x = [0._RKC, 0.5_RKC, 1._RKC], alpha = 2._RKC, beta = 3._RKC) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaCDF(x = [0._RKC, 0.5_RKC, 1._RKC], alpha = [0.1_RKC, 1._RKC, 10._RKC], beta = 3._RKC)")
    call disp%show( getBetaCDF(x = [0._RKC, 0.5_RKC, 1._RKC], alpha = [0.1_RKC, 1._RKC, 10._RKC], beta = 3._RKC) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 1000_IK
        real(RKC) :: X(NP)
        integer :: fileUnit, i

        call setLinSpace(X, 0._RKC, 1._RKC)
        open(newunit = fileUnit, file = "getBetaCDF.RK.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))" ) X(i), getBetaCDF(X(i), [.5_RKC, 5._RKC, .5_RKC, 5._RKC], [.5_RKC, 1.0_RKC, 5.0_RKC, 10._RKC])
        end do
        close(fileUnit)

    end block

end program example