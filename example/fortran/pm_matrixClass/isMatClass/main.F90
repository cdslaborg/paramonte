program example

    use pm_kind, only: SK, IK, RKG => RK
    use pm_io, only: display_type
    use pm_arrayChoice, only: getChoice
    use pm_distUnif, only: getUnifRand
    use pm_matrixClass, only: isMatClass, symmetric, hermitian, posdefmat
    use pm_matrixInit, only: getMatInit, uppLowDia
    use pm_distCov, only: getCovRand
    use pm_err, only: setAsserted
    use pm_val2str, only: getStr

    implicit none

    integer(IK) :: i, shape(2)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Check for Symmetric/Hermitian matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(2), allocatable :: choice(:), mat(:,:)
        do i = 1, 10
            call disp%skip
            call disp%show("choice = ['AA', 'BB'] ! example matrix element values.")
                            choice = ['AA', 'BB'] ! example matrix element values.
            call disp%show("mat = getMatInit(int(getUnifRand(2, 4, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getUnifRand(choice(1), choice(2)))")
                            mat = getMatInit(int(getUnifRand(2, 4, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getUnifRand(choice(1), choice(2)))
            call disp%show("mat")
            call disp%show( mat , deliml = """" )
            call disp%show("isMatClass(mat, class = symmetric)")
            call disp%show( isMatClass(mat, class = symmetric) )
            call disp%show("isMatClass(mat, class = hermitian)")
            call disp%show( isMatClass(mat, class = hermitian) )
        end do
    end block

    block
        integer, allocatable :: choice(:), mat(:,:)
        do i = 1, 10
            call disp%skip
            call disp%show("choice = [1, 2] ! example matrix element values.")
                            choice = [1, 2] ! example matrix element values.
            call disp%show("mat = getMatInit(int(getUnifRand(2, 4, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getUnifRand(choice(1), choice(2)))")
                            mat = getMatInit(int(getUnifRand(2, 4, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getUnifRand(choice(1), choice(2)))
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("isMatClass(mat, class = symmetric)")
            call disp%show( isMatClass(mat, class = symmetric) )
            call disp%show("isMatClass(mat, class = hermitian)")
            call disp%show( isMatClass(mat, class = hermitian) )
        end do
    end block

    block
        logical, allocatable :: choice(:), mat(:,:)
        do i = 1, 10
            call disp%skip
            call disp%show("choice = [.false., .true.] ! example matrix element values.")
                            choice = [.false., .true.] ! example matrix element values.
            call disp%show("mat = getMatInit(int(getUnifRand(2, 4, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getUnifRand(choice(1), choice(2)))")
                            mat = getMatInit(int(getUnifRand(2, 4, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getUnifRand(choice(1), choice(2)))
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("isMatClass(mat, class = symmetric)")
            call disp%show( isMatClass(mat, class = symmetric) )
            call disp%show("isMatClass(mat, class = hermitian)")
            call disp%show( isMatClass(mat, class = hermitian) )
        end do
    end block

    block
        complex, allocatable :: choice(:), mat(:,:)
        do i = 1, 10
            call disp%skip
            call disp%show("choice = [(1., -1.), (1., 0.), (1., 1.)] ! example matrix element values.")
                            choice = [(1., -1.), (1., 0.), (1., 1.)] ! example matrix element values.
            call disp%show("mat = getMatInit(int(spread(getUnifRand(2, 4), 1, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getChoice(choice))")
                            mat = getMatInit(int(spread(getUnifRand(2, 4), 1, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getChoice(choice))
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("isMatClass(mat, class = symmetric)")
            call disp%show( isMatClass(mat, class = symmetric) )
            call disp%show("isMatClass(mat, class = hermitian)")
            call disp%show( isMatClass(mat, class = hermitian) )
        end do
    end block

    block
        real, allocatable :: choice(:), mat(:,:)
        do i = 1, 10
            call disp%skip
            call disp%show("choice = [1., 2., 3.] ! example matrix element values.")
                            choice = [1., 2., 3.] ! example matrix element values.
            call disp%show("mat = getMatInit(int(spread(getUnifRand(2, 4), 1, 2), IK), uppLowDia, vupp = getUnifRand(choice(1), choice(2)), vlow = getUnifRand(choice(1), choice(2)), vdia = getUnifRand(choice(1), choice(2)))")
                            mat = getMatInit(int(spread(getUnifRand(2, 4), 1, 2), IK), uppLowDia, vupp = getUnifRand(choice(1), choice(2)), vlow = getUnifRand(choice(1), choice(2)), vdia = getUnifRand(choice(1), choice(2)))
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("isMatClass(mat, class = symmetric)")
            call disp%show( isMatClass(mat, class = symmetric) )
            call disp%show("isMatClass(mat, class = hermitian)")
            call disp%show( isMatClass(mat, class = hermitian) )
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Check for positive-definite matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        complex, allocatable :: choice(:), mat(:,:)
        do i = 1, 10
            call disp%skip
            call disp%show("choice = [(.5, -.5), (1., 0.), (.5, .5)] ! example matrix element values.")
                            choice = [(.5, -.5), (1., 0.), (.5, .5)] ! example matrix element values.
            call disp%show("mat = getMatInit(int(spread(getUnifRand(2, 4), 1, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getChoice(choice))")
                            mat = getMatInit(int(spread(getUnifRand(2, 4), 1, 2), IK), uppLowDia, vupp = getChoice(choice), vlow = getChoice(choice), vdia = getChoice(choice))
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("isMatClass(mat, class = hermitian)")
            call disp%show( isMatClass(mat, class = hermitian) )
            call disp%show("isMatClass(mat, class = posdefmat)")
            call disp%show( isMatClass(mat, class = posdefmat) )
        end do
    end block

    block
        use pm_kind, only: TKG => RKD
        real(TKG), allocatable :: choice(:), mat(:,:)
        do i = 1, 10
            call disp%skip
            call disp%show("mat = getCovRand(mold = 5._TKG, ndim = int(getUnifRand(2, 8), IK))")
                            mat = getCovRand(mold = 5._TKG, ndim = int(getUnifRand(2, 8), IK))
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("isMatClass(mat, class = symmetric)")
            call disp%show( isMatClass(mat, class = symmetric) )
            call disp%show("isMatClass(mat, class = hermitian)")
            call disp%show( isMatClass(mat, class = hermitian) )
            call disp%show("isMatClass(mat, class = posdefmat)")
            call disp%show( isMatClass(mat, class = posdefmat) )
            call disp%show("if (.not. isMatClass(mat, class = posdefmat)) error stop 'The output matrix from `getCovRand()` must be positive-definite. Please report this to developers.'")
                            if (.not. isMatClass(mat, class = posdefmat)) error stop 'The output matrix from `getCovRand()` must be positive-definite. Please report this to developers.'
        end do
    end block

end program example