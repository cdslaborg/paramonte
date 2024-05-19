program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_matrixChol, only: setChoLow
    use pm_distUnifEll, only: getUnifEllRand
    use pm_matrixClass, only: isMatClass, posdefmat
    use pm_distCov, only: getCovRand

    implicit none

    integer(IK) :: itry, ndim

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKG => RKS ! all real kinds are supported.
        real(TKG), allocatable :: rand(:,:)
        real(TKG) :: mold
        real(TKG) :: scale
        do itry = 1, 5
            call disp%skip()
            call disp%show("ndim = getUnifRand(2, 3)")
                            ndim = getUnifRand(2, 3)
            call disp%show("scale = getUnifRand(1, 10)")
                            scale = getUnifRand(1, 10)
            call disp%show("rand = getCovRand(mold, ndim)")
                            rand = getCovRand(mold, ndim)
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("rand = getCovRand(mold, ndim, scale)")
                            rand = getCovRand(mold, ndim, scale)
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%skip()
            call disp%skip()
            call disp%show("rand = getCovRand(mold = 0._TKG, scale = [(real(itry, TKG), itry = 1, 5)])")
                            rand = getCovRand(mold = 0._TKG, scale = [(real(itry, TKG), itry = 1, 5)])
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: TKG => CKS ! all real kinds are supported.
        complex(TKG), allocatable :: rand(:,:)
        complex(TKG) :: mold
        real(TKG) :: scale
        do itry = 1, 5
            call disp%skip()
            call disp%show("ndim = getUnifRand(2, 3)")
                            ndim = getUnifRand(2, 3)
            call disp%show("scale = getUnifRand(1, 10)")
                            scale = getUnifRand(1, 10)
            call disp%show("rand = getCovRand(mold, ndim)")
                            rand = getCovRand(mold, ndim)
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("rand = getCovRand(mold, ndim, scale)")
                            rand = getCovRand(mold, ndim, scale)
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%skip()
            call disp%skip()
            call disp%show("rand = getCovRand(mold = 0._TKG, scale = [(real(itry, TKG), itry = 1, 5)])")
                            rand = getCovRand(mold = 0._TKG, scale = [(real(itry, TKG), itry = 1, 5)])
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%skip()
        end do
    end block

end program example