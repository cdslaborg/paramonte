program example

    use pm_kind, only: SK, IK
    use pm_cosmology, only: getDisComTransNormedWU10
    use pm_io, only: display_type

    implicit none

    type(display_type)      :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Transverse Comoving Distance in units of Hubble Distance.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getDisComTransNormedWU10(zplus1 = 1.)")
    call disp%show( getDisComTransNormedWU10(zplus1 = 1.) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDisComTransNormedWU10(zplus1 = 1.1)")
    call disp%show( getDisComTransNormedWU10(zplus1 = 1.1) )
    call disp%skip()

    call disp%skip()
    call disp%show("getDisComTransNormedWU10(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.6)")
    call disp%show( getDisComTransNormedWU10(zplus1 = [ 2., 3.], omegaM = 0.4, omegaL = 0.6) )
    call disp%skip()

    ! Generate both the cosmic rate and the rate density.

    block

        use pm_val2str, only: getStr
        use pm_arraySpace, only: getLinSpace
        use pm_arraySpace, only: getLogSpace

        real, allocatable   :: zplus1(:), OmegaM(:), OmegaL(:)
        integer :: fileUnit, i

        OmegaM = [0.2, 0.3, 0.4]
        OmegaL = 1. - OmegaM
        zplus1 = 1. + getLogSpace(log(0.0001), log(10000.), 500_IK)

        open(newunit = fileUnit, file = "getDisComTransNormedWU10.RK.txt")
        write(fileUnit, "(*(g0,:,','))") "z", ("DisLum_"//getStr(OmegaM(i),"(g0.1)")//"_"//getStr(OmegaL(i),"(g0.1)"), i = 1, size(OmegaM))
        do i = 1, size(zplus1)
            write(fileUnit, "(*(g0,:,','))") zplus1(i) - 1., getDisComTransNormedWU10(zplus1(i), OmegaM, OmegaL)
        end do
        close(fileUnit)

    end block

end program example