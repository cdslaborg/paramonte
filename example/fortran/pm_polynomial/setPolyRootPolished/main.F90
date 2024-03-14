program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_polynomial, only: getPolyRoot
    use pm_polynomial, only: setPolyRootPolished
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: niter, iroot
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => RKH ! all processor kinds are supported.
        real(TKC), parameter :: noise = sqrt(epsilon(0._TKC))
        complex(TKC), allocatable :: root(:)
        real(TKC), allocatable :: coef(:)
        call disp%skip()
        call disp%show("coef = [real(TKC) :: -4, 1, -4, 1]")
                        coef = [real(TKC) :: -4, 1, -4, 1]
        call disp%show("root = getPolyRoot(coef)")
                        root = getPolyRoot(coef)
        call disp%show("if (size(root) == 0_IK) error stop 'no roots could be found.'")
                        if (size(root) == 0_IK) error stop 'no roots could be found.'
        call disp%show("root")
        call disp%show( root )
        call disp%show("root = root + (noise, noise)")
                        root = root + (noise, noise)
        call disp%show("root")
        call disp%show( root )
        do iroot = 1, size(root)
            call disp%show("iroot")
            call disp%show( iroot )
            call disp%show("root(iroot)")
            call disp%show( root(iroot) )
            call disp%show("call setPolyRootPolished(root(iroot), niter, coef)")
                            call setPolyRootPolished(root(iroot), niter, coef)
            call disp%show("niter")
            call disp%show( niter )
            call disp%show("if (0 < niter) call disp%show(root(iroot)) ! polished")
                            if (0 < niter) call disp%show(root(iroot)) ! polished
        end do
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKH ! all processor kinds are supported.
        real(TKC), parameter :: noise = sqrt(epsilon(0._TKC))
        complex(TKC), allocatable :: coef(:), root(:)
        call disp%skip()
        call disp%show("coef = cmplx([1, 3, 0, 2], -[1, 3, 0, 2], TKC)")
                        coef = cmplx([1, 3, 0, 2], -[1, 3, 0, 2], TKC)
        call disp%show("root = getPolyRoot(coef)")
                        root = getPolyRoot(coef)
        call disp%show("if (size(root) == 0_IK) error stop 'no roots could be found.'")
                        if (size(root) == 0_IK) error stop 'no roots could be found.'
        call disp%show("root")
        call disp%show( root )
        call disp%show("root = root + (noise, noise)")
                        root = root + (noise, noise)
        call disp%show("root")
        call disp%show( root )
        do iroot = 1, size(root)
        call disp%show("iroot")
        call disp%show( iroot )
        call disp%show("root(iroot)")
        call disp%show( root(iroot) )
        call disp%show("call setPolyRootPolished(root(iroot), niter, coef)")
                        call setPolyRootPolished(root(iroot), niter, coef)
        call disp%show("niter")
        call disp%show( niter )
        call disp%show("if (0 < niter) call disp%show(root(iroot)) ! polished")
                        if (0 < niter) call disp%show(root(iroot)) ! polished
        end do
        call disp%skip()
    end block

end program example