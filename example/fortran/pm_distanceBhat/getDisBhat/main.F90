program example

    use pm_kind, only: SK, IK
    use pm_arrayRange, only: getRange
    use pm_arrayResize, only: setResized
    use pm_distanceBhat, only: getDisBhat
    use pm_distGeomCyclic, only: getGeomCyclicLogPMF
    use pm_distNorm, only: getNormLogPDF
    use pm_distCov, only: getCovRand
    use pm_io, only: display_type

    implicit none

    integer(IK) :: ndim, npnt, nsam, isam

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Bhattacharyya distance squared for real-valued PMFs.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: RKG => RKS ! all other real kinds are also supported.
        real(RKG), allocatable :: p(:), q(:), stepSuccess(:)
        integer(IK) :: period
        real(RKG) :: bhat

        call disp%skip()
        call disp%show("period = 10;")
                        period = 10;
        call disp%show("stepSuccess = getRange(1, period)")
                        stepSuccess = getRange(1, period)
        call disp%show("stepSuccess")
        call disp%show( stepSuccess )
        call disp%show("p = exp(getGeomCyclicLogPMF(stepSuccess, probSuccess = .1_RKG, period = period))")
                        p = exp(getGeomCyclicLogPMF(stepSuccess, probSuccess = .1_RKG, period = period))
        call disp%show("p")
        call disp%show( p )
        call disp%show("q = p")
                        q = p
        call disp%show("q")
        call disp%show( q )
        call disp%show("bhat = getDisBhat(p, q)")
                        bhat = getDisBhat(p, q)
        call disp%show("bhat")
        call disp%show( bhat )
        call disp%skip()

        call disp%skip()
        call disp%show("period = 10;")
                        period = 10;
        call disp%show("stepSuccess = getRange(1, period)")
                        stepSuccess = getRange(1, period)
        call disp%show("stepSuccess")
        call disp%show( stepSuccess )
        call disp%show("p = exp(getGeomCyclicLogPMF(stepSuccess, probSuccess = .1_RKG, period = period))")
                        p = exp(getGeomCyclicLogPMF(stepSuccess, probSuccess = .1_RKG, period = period))
        call disp%show("p")
        call disp%show( p )
        call disp%show("q = exp(getGeomCyclicLogPMF(stepSuccess, probSuccess = .9_RKG, period = period))")
                        q = exp(getGeomCyclicLogPMF(stepSuccess, probSuccess = .9_RKG, period = period))
        call disp%show("q")
        call disp%show( q )
        call disp%show("bhat = getDisBhat(p, q)")
                        bhat = getDisBhat(p, q)
        call disp%show("bhat")
        call disp%show( bhat )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Bhattacharyya distance squared for real-valued PDFs.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: RKG => RKD ! all other real kinds are also supported.
        real(RKG) :: bhat

        call disp%skip()
        call disp%show("bhat = getDisBhat(getp, getq)")
                        bhat = getDisBhat(getp, getq)
        call disp%show("bhat")
        call disp%show( bhat )
        call disp%skip()
        call disp%show("bhat = getDisBhat(getp, getq, lb = -huge(0._RKG), ub = +huge(0._RKG))")
                        bhat = getDisBhat(getp, getq, lb = -huge(0._RKG), ub = +huge(0._RKG))
        call disp%show("bhat")
        call disp%show( bhat )
        call disp%skip()

    end block

contains

    function getp(x) result(pdf)
        use pm_kind, only: RKG => RKD ! all other real kinds are also supported.
        real(RKG), intent(in) :: x
        real(RKG) :: pdf
        pdf = exp(getNormLogPDF(x))
    end function

    function getq(x) result(pdf)
        use pm_kind, only: RKG => RKD ! all other real kinds are also supported.
        real(RKG), intent(in) :: x
        real(RKG) :: pdf
        pdf = exp(getNormLogPDF(x, mu = 1._RKG))
    end function

end program example