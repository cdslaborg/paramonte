program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_container, only: csp_type
    use pm_distUnif, only: getUnifRand
    use pm_matrixClass, only: isMatClass, posdefmat
    use pm_distChol, only: rngf_type, xoshiro256ssw_type
    use pm_distChol, only: setCholRand, uppDia, lowDia, subset_type
    use pm_arrayRefill, only: setRefilled

    implicit none

    integer(IK) :: isub, itry, ndim

    type(csp_type) :: subset(2)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    subset = [csp_type(uppDia), csp_type(lowDia)]
    do isub = 1, size(subset)

        block
            use pm_kind, only: TKG => RKS ! all real kinds are supported.
            real(TKG), allocatable :: rand(:,:), cov(:,:)
            real(TKG) :: mold
            do itry = 1, 5
                call disp%skip()
                call disp%show("ndim = getUnifRand(2, 5); call setRefilled(rand, 0._TKG, [ndim, ndim])")
                                ndim = getUnifRand(2, 5); call setRefilled(rand, 0._TKG, [ndim, ndim])
                call disp%show("[same_type_as(subset(isub)%val, uppDia), same_type_as(subset(isub)%val, lowDia)]")
                call disp%show( [same_type_as(subset(isub)%val, uppDia), same_type_as(subset(isub)%val, lowDia)] )
                call disp%show("call setCholRand(rngf_type(), rand, subset(isub)%val)")
                                call setCholRand(rngf_type(), rand, subset(isub)%val)
                call disp%show("rand")
                call disp%show( rand )
                call disp%show("cov = matmul(rand, transpose(rand))")
                                cov = matmul(rand, transpose(rand))
                call disp%show("isMatClass(cov, posdefmat)")
                call disp%show( isMatClass(cov, posdefmat) )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKG => CKS ! all complex kinds are supported.
            complex(TKG), allocatable :: rand(:,:), cov(:,:)
            complex(TKG) :: mold
            do itry = 1, 5
                call disp%skip()
                call disp%show("ndim = getUnifRand(2, 5); call setRefilled(rand, cmplx(0, 0, TKG), [ndim, ndim])")
                                ndim = getUnifRand(2, 5); call setRefilled(rand, cmplx(0, 0, TKG), [ndim, ndim])
                call disp%show("[same_type_as(subset(isub)%val, uppDia), same_type_as(subset(isub)%val, lowDia)]")
                call disp%show( [same_type_as(subset(isub)%val, uppDia), same_type_as(subset(isub)%val, lowDia)] )
                call disp%show("call setCholRand(rngf_type(), rand, subset(isub)%val)")
                                call setCholRand(rngf_type(), rand, subset(isub)%val)
                call disp%show("rand")
                call disp%show( rand )
                call disp%show("cov = matmul(rand, conjg(transpose(rand)))")
                                cov = matmul(rand, conjg(transpose(rand)))
                call disp%show("isMatClass(cov, posdefmat)")
                call disp%show( isMatClass(cov, posdefmat) )
                call disp%skip()
            end do
        end block

    end do

end program example