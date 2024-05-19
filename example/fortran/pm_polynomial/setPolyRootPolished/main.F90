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
        use pm_kind, only: TKG => RKH ! all processor kinds are supported.
        real(TKG), parameter :: noise = sqrt(epsilon(0._TKG))
        complex(TKG), allocatable :: root(:)
        real(TKG), allocatable :: coef(:)
        call disp%skip()
        call disp%show("coef = [real(TKG) :: -4, 1, -4, 1]")
                        coef = [real(TKG) :: -4, 1, -4, 1]
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
        use pm_kind, only: TKG => RKH ! all processor kinds are supported.
        real(TKG), parameter :: noise = sqrt(epsilon(0._TKG))
        complex(TKG), allocatable :: coef(:), root(:)
        call disp%skip()
        call disp%show("coef = cmplx([1, 3, 0, 2], -[1, 3, 0, 2], TKG)")
                        coef = cmplx([1, 3, 0, 2], -[1, 3, 0, 2], TKG)
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