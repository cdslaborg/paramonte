program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK
    use pm_io, only: field_type
    use pm_io, only: display_type
    use pm_matrixChol, only: setChoLow
    use pm_distUnif, only: getUnifRand, rngf
    use pm_distUnifEll, only: getUnifEllRand
    use pm_distCov, only: setCovRand, dvine, onion
    use pm_matrixClass, only: isMatClass, posdefmat, hermitian
    use pm_arrayResize, only: setResized
    use pm_matrixDet, only: getMatDet

    implicit none

    integer(IK) :: itry, ndim

    type(display_type) :: disp
    disp = display_type(file = SK_"main.out.F90", format = field_type(complex = SK_"math"))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Gram method for real covariance.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKC => RKS ! all kinds are supported.
        real(TKC), allocatable :: scale(:)
        real(TKC), allocatable :: rand(:,:)
        do itry = 1, 5

            call disp%skip()
            call disp%show("ndim = getUnifRand(3, 9)")
                            ndim = getUnifRand(3, 9)
            call disp%show("ndim")
            call disp%show( ndim )
            call disp%show("call setResized(rand, [ndim, ndim])")
                            call setResized(rand, [ndim, ndim])
            call disp%show("scale = getUnifRand(1, 10, ndim)")
                            scale = getUnifRand(1, 10, ndim)
            call disp%show("scale")
            call disp%show( scale )
            call disp%skip()

            call disp%show("call setCovRand(rngf, rand)")
                            call setCovRand(rngf, rand)
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )

            call disp%show("call setCovRand(rngf, rand, scale(1))")
                            call setCovRand(rngf, rand, scale(1))
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )
            call disp%skip()

            call disp%show("call setCovRand(rngf, rand, scale)")
                            call setCovRand(rngf, rand, scale)
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )
            call disp%skip()

        end do

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Gram method for complex covariance.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKC => RKD ! all kinds are supported.
        real(TKC), allocatable :: scale(:)
        complex(TKC), allocatable :: rand(:,:)
        do itry = 1, 5

            call disp%skip()
            call disp%show("ndim = getUnifRand(3, 9)")
                            ndim = getUnifRand(3, 9)
            call disp%show("ndim")
            call disp%show( ndim )
            call disp%show("call setResized(rand, [ndim, ndim])")
                            call setResized(rand, [ndim, ndim])
            call disp%show("scale = getUnifRand(1, 10, ndim)")
                            scale = getUnifRand(1, 10, ndim)
            call disp%show("scale")
            call disp%show( scale )
            call disp%skip()

            call disp%show("call setCovRand(rngf, rand)")
                            call setCovRand(rngf, rand)
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )

            call disp%show("call setCovRand(rngf, rand, scale(1))")
                            call setCovRand(rngf, rand, scale(1))
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )
            call disp%skip()

            call disp%show("call setCovRand(rngf, rand, scale)")
                            call setCovRand(rngf, rand, scale)
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )
            call disp%skip()

        end do

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Dvine and Onion methods.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKC => RKS ! all kinds are supported.
        real(TKC), allocatable :: rand(:,:)
        real(TKC)   :: eta, scale
        do itry = 1, 10

            call disp%skip()
            call disp%show("eta = getUnifRand(1, 10)")
                            eta = getUnifRand(1, 10)
            call disp%show("ndim = getUnifRand(2, 5)")
                            ndim = getUnifRand(2, 5)
            call disp%show("call setResized(rand, [ndim, ndim])")
                            call setResized(rand, [ndim, ndim])
            call disp%show("scale = getUnifRand(1, 10)")
                            scale = getUnifRand(1, 10)
            call disp%skip()

            call disp%show("call setCovRand(rngf, rand, dvine, eta)")
                            call setCovRand(rngf, rand, dvine, eta)
            call disp%show("onion%info")
            call disp%show( onion%info )
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )

            call disp%show("call setCovRand(rngf, rand, onion, eta)")
                            call setCovRand(rngf, rand, onion, eta)
            call disp%show("onion%info")
            call disp%show( onion%info )
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )

            call disp%show("call setCovRand(rngf, rand, dvine, eta, scale)")
                            call setCovRand(rngf, rand, dvine, eta, scale)
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )
            call disp%skip()

            call disp%show("call setCovRand(rngf, rand, onion, eta, scale)")
                            call setCovRand(rngf, rand, onion, eta, scale)
            call disp%show("onion%info")
            call disp%show( onion%info )
            call disp%show("rand")
            call disp%show( rand )
            call disp%show("isMatClass(rand, posdefmat)")
            call disp%show( isMatClass(rand, posdefmat) )
            call disp%show("isMatClass(rand, hermitian)")
            call disp%show( isMatClass(rand, hermitian) )
            call disp%show("getMatDet(rand)")
            call disp%show( getMatDet(rand) )
            call disp%skip()

        end do

        call disp%skip()
        call disp%show("ndim = getUnifRand(2, 10)")
                        ndim = getUnifRand(2, 10)
        call disp%show("call setResized(rand, [ndim, ndim])")
                        call setResized(rand, [ndim, ndim])
        call disp%skip()

        call disp%show("call setCovRand(rngf, rand, dvine, eta = 0._TKC, scale = [(real(itry, TKC), itry = 1, ndim)])")
                        call setCovRand(rngf, rand, dvine, eta = 0._TKC, scale = [(real(itry, TKC), itry = 1, ndim)])
        call disp%show("rand")
        call disp%show( rand )
        call disp%show("isMatClass(rand, posdefmat)")
        call disp%show( isMatClass(rand, posdefmat) )
        call disp%skip()

        call disp%show("call setCovRand(rngf, rand, onion, eta = 0._TKC, scale = [(real(itry, TKC), itry = 1, ndim)])")
                        call setCovRand(rngf, rand, onion, eta = 0._TKC, scale = [(real(itry, TKC), itry = 1, ndim)])
        call disp%show("onion%info")
        call disp%show( onion%info )
        call disp%show("rand")
        call disp%show( rand )
        call disp%show("isMatClass(rand, posdefmat)")
        call disp%show( isMatClass(rand, posdefmat) )
        call disp%skip()

    end block

end program example