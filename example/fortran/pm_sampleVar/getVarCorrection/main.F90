program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_sampleVar, only: getVarCorrection

    implicit none

    integer(IK) :: itry, ntry = 10
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: RKC => RKS ! All other real types are also supported.
        integer(RKC), allocatable :: weight(:)
        real(RKC) :: correction
        do itry = 1, ntry
            call disp%skip()
            call disp%show("weight = getUnifRand(1, 9, getUnifRand(2_IK, 20_IK))")
                            weight = getUnifRand(1, 9, getUnifRand(2_IK, 20_IK))
            call disp%show("weight")
            call disp%show( weight )
            call disp%show("correction = getVarCorrection(real(sum(weight), RKC))")
                            correction = getVarCorrection(real(sum(weight), RKC))
            call disp%show("correction")
            call disp%show( correction )
            call disp%show("correction = getVarCorrection(real(sum(weight), RKC), real(sum(weight**2), RKC))")
                            correction = getVarCorrection(real(sum(weight), RKC), real(sum(weight**2), RKC))
            call disp%show("correction")
            call disp%show( correction )
            call disp%skip()
        end do
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_kind, only: RKC => RK
        use pm_arraySpace, only: getLogSpace
        use pm_io, only: getErrTableWrite, trans
        real(RKC), allocatable :: sampleSize(:)
        sampleSize = getLogSpace(log(1.1_RKC), log(1001._RKC), 1000_IK)
        if (0 /= getErrTableWrite("getVarCorrection.RK.txt", reshape([sampleSize, getVarCorrection(sampleSize)], [size(sampleSize), 2]), header = SK_"Sample Size,Bessel Correction")) error stop 'table write failed.'
    end block

end program example