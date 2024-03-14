program example

    use pm_kind, only: SK, IK
    use pm_cosmology, only: getVolComDiffNormed
    use pm_io, only: display_type

    implicit none

    type(display_type)      :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Comoving Volume Element in units of Hubble Volume.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getVolComDiffNormed(zplus1 = 1., reltol = sqrt(epsilon(0.)))")
    call disp%show( getVolComDiffNormed(zplus1 = 1., reltol = sqrt(epsilon(0.))) )
    call disp%skip()

    call disp%skip()
    call disp%show("getVolComDiffNormed(zplus1 = 1.1, reltol = sqrt(epsilon(0.)))")
    call disp%show( getVolComDiffNormed(zplus1 = 1.1, reltol = sqrt(epsilon(0.))) )
    call disp%skip()

    call disp%skip()
    call disp%show("getVolComDiffNormed(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.6, reltol = 0.0001)")
    call disp%show( getVolComDiffNormed(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.6, reltol = 0.0001) )
    call disp%skip()

    call disp%skip()
    call disp%show("getVolComDiffNormed(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.4, omegaR = 0.2, reltol = 0.0001)")
    call disp%show( getVolComDiffNormed(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.4, omegaR = 0.2, reltol = 0.0001) )
    call disp%skip()

    call disp%skip()
    call disp%show("getVolComDiffNormed(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.4, omegaR = 0.4, omegaK = -0.2, sqrtAbsOmegaK = sqrt(0.2), reltol = 0.0001)")
    call disp%show( getVolComDiffNormed(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.4, omegaR = 0.4, omegaK = -0.2, sqrtAbsOmegaK = sqrt(0.2), reltol = 0.0001) )
    call disp%skip()

    ! Generate both the cosmic rate and the rate density.

    block

        use pm_val2str, only: getStr
        use pm_arraySpace, only: getLinSpace
        use pm_arraySpace, only: getLogSpace

        real, allocatable   :: zplus1(:), OmegaM(:), OmegaL(:)
        integer :: fileUnit, i

        OmegaM = [1., 0.3, .05]
        OmegaL = 1. - OmegaM
        zplus1 = 1. + getLogSpace(log(0.0001), log(10000.), 500_IK)

        open(newunit = fileUnit, file = "getVolComDiffNormed.RK.txt")
        write(fileUnit, "(*(g0,:,','))") "z", ("DisLum_"//getStr(OmegaM(i),"(g0.1)")//"_"//getStr(OmegaL(i),"(g0.1)"), i = 1, size(OmegaM))
        do i = 1, size(zplus1)
            write(fileUnit, "(*(g0,:,','))") zplus1(i) - 1., getVolComDiffNormed(zplus1(i), OmegaM, OmegaL, sqrt(epsilon(0.)))
        end do
        close(fileUnit)

    end block

end program example