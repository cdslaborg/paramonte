program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBeta, only: getBetaCDF

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getBetaCDF(x = 0._RKG, alpha = 2._RKG, beta = 3._RKG)")
    call disp%show( getBetaCDF(x = 0._RKG, alpha = 2._RKG, beta = 3._RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaCDF(x = .5_RKG, alpha = 2._RKG, beta = 3._RKG)")
    call disp%show( getBetaCDF(x = .5_RKG, alpha = 2._RKG, beta = 3._RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaCDF(x = 1._RKG, alpha = 2._RKG, beta = 3._RKG)")
    call disp%show( getBetaCDF(x = 1._RKG, alpha = 2._RKG, beta = 3._RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaCDF(x = [0._RKG, 0.5_RKG, 1._RKG], alpha = 2._RKG, beta = 3._RKG)")
    call disp%show( getBetaCDF(x = [0._RKG, 0.5_RKG, 1._RKG], alpha = 2._RKG, beta = 3._RKG) )
    call disp%skip()

    call disp%skip()
    call disp%show("getBetaCDF(x = [0._RKG, 0.5_RKG, 1._RKG], alpha = [0.1_RKG, 1._RKG, 10._RKG], beta = 3._RKG)")
    call disp%show( getBetaCDF(x = [0._RKG, 0.5_RKG, 1._RKG], alpha = [0.1_RKG, 1._RKG, 10._RKG], beta = 3._RKG) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 1000_IK
        real(RKG) :: X(NP)
        integer :: fileUnit, i

        call setLinSpace(X, 0._RKG, 1._RKG)
        open(newunit = fileUnit, file = "getBetaCDF.RK.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))") X(i), getBetaCDF(X(i), [.5_RKG, 5._RKG, .5_RKG, 5._RKG], [.5_RKG, 1.0_RKG, 5.0_RKG, 10._RKG])
        end do
        close(fileUnit)

    end block

end program example