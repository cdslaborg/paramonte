program example

    use pm_kind, only: SK, IK
    use pm_cosmology, only: getDisComTransNormed
    use pm_cosmology, only: getHubbleParamNormedSq
    use pm_cosmology, only: setVolComDiffNormed
    use pm_io, only: display_type

    implicit none

    real :: VolComDiffNormed(3)
    type(display_type)      :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Comoving Volume Element in units of Hubble Volume.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setVolComDiffNormed(VolComDiffNormed(1), getDisComTransNormed(zplus1 = 1., reltol = sqrt(epsilon(0.)))**2, sqrt(getHubbleParamNormedSq(zplus1 = 1.)))")
                    call setVolComDiffNormed(VolComDiffNormed(1), getDisComTransNormed(zplus1 = 1., reltol = sqrt(epsilon(0.)))**2, sqrt(getHubbleParamNormedSq(zplus1 = 1.)))
    call disp%show("VolComDiffNormed(1)")
    call disp%show( VolComDiffNormed(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setVolComDiffNormed(VolComDiffNormed(1), getDisComTransNormed(zplus1 = 1.1, reltol = sqrt(epsilon(0.)))**2, sqrt(getHubbleParamNormedSq(zplus1 = 1.1)))")
                    call setVolComDiffNormed(VolComDiffNormed(1), getDisComTransNormed(zplus1 = 1.1, reltol = sqrt(epsilon(0.)))**2, sqrt(getHubbleParamNormedSq(zplus1 = 1.1)))
    call disp%show("VolComDiffNormed(1)")
    call disp%show( VolComDiffNormed(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setVolComDiffNormed(VolComDiffNormed(1:2), getDisComTransNormed(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.6, reltol = 0.0001)**2, sqrt(getHubbleParamNormedSq(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.6)))")
                    call setVolComDiffNormed(VolComDiffNormed(1:2), getDisComTransNormed(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.6, reltol = 0.0001)**2, sqrt(getHubbleParamNormedSq(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.6)))
    call disp%show("VolComDiffNormed(1:2)")
    call disp%show( VolComDiffNormed(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setVolComDiffNormed(VolComDiffNormed(1:2), getDisComTransNormed(zplus1 = [ 2., 3.], omegaM = 0.2, omegaL = 0.6, omegaR = 0.2, reltol = 0.0001)**2, sqrt(getHubbleParamNormedSq(zplus1 = [ 2., 3.], omegaM = 0.2, omegaL = 0.6, omegaR = 0.2)))")
                    call setVolComDiffNormed(VolComDiffNormed(1:2), getDisComTransNormed(zplus1 = [ 2., 3.], omegaM = 0.2, omegaL = 0.6, omegaR = 0.2, reltol = 0.0001)**2, sqrt(getHubbleParamNormedSq(zplus1 = [ 2., 3.], omegaM = 0.2, omegaL = 0.6, omegaR = 0.2)))
    call disp%show("VolComDiffNormed(1:2)")
    call disp%show( VolComDiffNormed(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setVolComDiffNormed(VolComDiffNormed(1:2), getDisComTransNormed(zplus1 = [ 2., 3.], omegaM = 0.2, omegaL = 0.6, omegaR = 0.4, omegaK = -0.2, sqrtAbsOmegaK = sqrt(0.2), reltol = 0.0001)**2, sqrt(getHubbleParamNormedSq(zplus1 = [ 2., 3.], omegaM = 0.2, omegaL = 0.6, omegaR = 0.4, omegaK = -0.2)))")
                    call setVolComDiffNormed(VolComDiffNormed(1:2), getDisComTransNormed(zplus1 = [ 2., 3.], omegaM = 0.2, omegaL = 0.6, omegaR = 0.4, omegaK = -0.2, sqrtAbsOmegaK = sqrt(0.2), reltol = 0.0001)**2, sqrt(getHubbleParamNormedSq(zplus1 = [ 2., 3.], omegaM = 0.2, omegaL = 0.6, omegaR = 0.4, omegaK = -0.2)))
    call disp%show("VolComDiffNormed(1:2)")
    call disp%show( VolComDiffNormed(1:2) )
    call disp%skip()

    ! Generate both the cosmic rate and the rate density.

    block

        use pm_val2str, only: getStr
        use pm_arraySpace, only: getLinSpace
        use pm_arraySpace, only: getLogSpace

        real, allocatable   :: zplus1(:), OmegaM(:), OmegaL(:), VolComDiffNormed(:)
        integer :: fileUnit, i

        OmegaM = [1., 0.3, .05]
        OmegaL = 1. - OmegaM
        allocate(VolComDiffNormed, mold = OmegaM)
        zplus1 = 1. + getLogSpace(log(0.0001), log(10000.), 500_IK)

        open(newunit = fileUnit, file = "setVolComDiffNormed.RK.txt")
        write(fileUnit, "(*(g0,:,','))") "z", ("DisLum_"//getStr(OmegaM(i),"(g0.1)")//"_"//getStr(OmegaL(i),"(g0.1)"), i = 1, size(OmegaM))
        do i = 1, size(zplus1)
            call setVolComDiffNormed(VolComDiffNormed, getDisComTransNormed(zplus1(i), OmegaM, OmegaL, 0.0001)**2, sqrt(getHubbleParamNormedSq(zplus1(i), OmegaM, OmegaL)))
            write(fileUnit, "(*(g0,:,','))") zplus1(i) - 1., VolComDiffNormed       
        end do
        close(fileUnit)

    end block

end program example