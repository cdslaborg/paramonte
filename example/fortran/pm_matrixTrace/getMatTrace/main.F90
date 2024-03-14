program example

    use iso_fortran_env, only: output_unit
    use pm_kind, only: SK, IK, LK, CK, RK
    use pm_matrixSubset, only: dia, uppDia, lowDia, uppLow, uppLowDia
    use pm_matrixCopy, only: getMatCopy, rdpack, rfpack, lfpack
    use pm_matrixTrace, only: getMatTrace
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    integer(IK) :: ndim, itry, ntry = 5

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    block
        integer, allocatable :: mat(:,:), rfp(:,:), lfp(:) ! all integer kinds are supported.
        do itry = 1, ntry
            call disp%skip()
            call disp%show("ndim = getUnifRand(1, 7)")
                            ndim = getUnifRand(1, 7)
            call disp%show("mat = getUnifRand(1, 9, ndim, ndim)")
                            mat = getUnifRand(1, 9, ndim, ndim)
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("getMatTrace(mat)")
            call disp%show( getMatTrace(mat) )
            call disp%show("getMatTrace(getMatCopy(rdpack, mat, rdpack, subset = uppDia))")
            call disp%show( getMatTrace(getMatCopy(rdpack, mat, rdpack, subset = uppDia)) )
            call disp%show("getMatTrace(getMatCopy(rdpack, mat, rdpack, subset = lowDia))")
            call disp%show( getMatTrace(getMatCopy(rdpack, mat, rdpack, subset = lowDia)) )
            call disp%show("getMatTrace(getMatCopy(lfpack, mat, rdpack, lowDia), lfpack, lowDia)")
            call disp%show( getMatTrace(getMatCopy(lfpack, mat, rdpack, lowDia), lfpack, lowDia) )
            call disp%skip()
        end do
    end block

end program example