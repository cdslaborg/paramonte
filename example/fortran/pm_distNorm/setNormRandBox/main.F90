program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RK ! all real kinds are supported.
    use pm_distNorm, only: setNormRandBox
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    real(RK), dimension(NP) :: rand

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    !call setLinSpace(mean, x1 = -5._RK, x2 = +5._RK)
    !call setLogSpace(std, logx1 = log(0.1_RK), logx2 = log(10._RK))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random numbers from the (Standard) Normal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Normal random number from a Standard Normal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call random_number(rand(1:2))")
                    call random_number(rand(1:2))
    call disp%show("call setNormRandBox(rand(1), rand(2))")
                    call setNormRandBox(rand(1), rand(2))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Box-Muller Basic (Trigonometric) method.

    block
        integer(IK) :: fileUnit, i
        integer(IK), parameter :: NP = 5000_IK
        real(RK), dimension(NP) :: rand1, rand2, rand3
        call random_number(rand1)
        call random_number(rand2)
        call random_number(rand3)
        call setNormRandBox(rand1(1::2), rand1(2::2)); rand1 = rand1 * 3.0_RK + 2._RK
        call setNormRandBox(rand2(1::2), rand2(2::2)); rand2 = rand2 * 1.0_RK + 0._RK
        call setNormRandBox(rand3(1::2), rand3(2::2)); rand3 = rand3 * 1.0_RK - 5._RK
        open(newunit = fileUnit, file = "setNormRandBox.RK.txt")
        write(fileUnit,"(3(g0,:,' '))") ( rand1(i) &
                                        , rand2(i) &
                                        , rand3(i) &
                                        , i = 1,NP &
                                        )
        close(fileUnit)
    end block

    ! Box-Muller Polar (Rejection) method.

    block
        logical :: failed
        integer(IK) :: fileUnit, i, j
        integer(IK), parameter :: NP = 5000_IK
        real(RK), dimension(2) :: rand1, rand2, rand3
        open(newunit = fileUnit, file = "setNormRandBox.RK.txt")
        do i = 1, NP
            do
                call random_number(rand1(1:2))
                call setNormRandBox(rand1(1), rand1(2), failed = failed)
                if (.not. failed) exit
            end do
            rand1(1:2) = rand1(1:2) * 3 + 2
            do
                call random_number(rand2(1:2))
                call setNormRandBox(rand2(1), rand2(2), failed = failed)
                if (.not. failed) exit
            end do
            do
                call random_number(rand3(1:2))
                call setNormRandBox(rand3(1), rand3(2), failed = failed)
                if (.not. failed) exit
            end do
            rand3(1:2) = rand3(1:2) - 5._RK
            do j = 1, 2
                write(fileUnit,"(3(g0,:,' '))") rand1(j), rand2(j), rand3(j)
            end do
        end do
        close(fileUnit)
    end block

end program example