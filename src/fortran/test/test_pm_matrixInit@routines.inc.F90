!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This file contains test implementations of the procedure of [pm_matrixInit](@ref pm_matrixInit).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMatInitULD_D2_LK_ENABLED || setMatInitULD_D2_LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     !(getMatInitULD_D2_ENABLED || setMatInitULD_D2_ENABLED)
#error  "Unrecognized Interface."
#endif

        use pm_distUnif, only: getUnifRand, setUnifRand
        use pm_arrayReverse, only: setReversed
        use pm_val2str, only: getStr
        use pm_kind, only: IK, LK

        integer(IK) :: itest, MatShape(2), nrow, ncol, roff, coff, doff, maxdim
        logical(LK) :: isScalarDia

#if     getMatInitULD_D2_SK_ENABLED || setMatInitULD_D2_SK_ENABLED
        character(2,SKG), allocatable   :: Mat_ref(:,:), mat(:,:), vlow, vupp, vdia(:)
        vlow = SKG_"aa"; vupp = SKG_"zz"
#elif   getMatInitULD_D2_IK_ENABLED || setMatInitULD_D2_IK_ENABLED
        integer(IKG)    , allocatable   :: Mat_ref(:,:), mat(:,:), vlow, vupp, vdia(:)
        vlow = -1_IKG; vupp = +1_IKG
#elif   getMatInitULD_D2_LK_ENABLED || setMatInitULD_D2_LK_ENABLED
        logical(LKG)    , allocatable   :: Mat_ref(:,:), mat(:,:), vlow, vupp, vdia(:)
        vlow = .false._LKG; vupp = .true._LKG
#elif   getMatInitULD_D2_CK_ENABLED || setMatInitULD_D2_CK_ENABLED
        complex(CKG)    , allocatable   :: Mat_ref(:,:), mat(:,:), vlow, vupp, vdia(:)
        vlow = (-1._CKG, -1._CKG); vupp = (1._CKG, 1._CKG)
#elif   getMatInitULD_D2_RK_ENABLED || setMatInitULD_D2_RK_ENABLED
        real(RKG)       , allocatable   :: Mat_ref(:,:), mat(:,:), vlow, vupp, vdia(:)
        vlow = -1._RKG; vupp = 1._RKG
#else
#error  "Unrecognized Interface."
#endif

        assertion = .true._LK
        maxdim = 10_IK

        do itest = 1, 200

            call setUnifRand(MatShape, 1_IK, maxdim)
            nrow = getUnifRand(1_IK, MatShape(1))
            ncol = getUnifRand(1_IK, MatShape(2))
            roff = getUnifRand(0_IK, MatShape(1) - nrow)
            coff = getUnifRand(0_IK, MatShape(2) - ncol)
            doff = getUnifRand(min(0_IK, 1_IK - nrow), max(0_IK, ncol - 1_IK))
            if (allocated(Mat_ref)) deallocate(Mat_ref)
            allocate(Mat_ref(MatShape(1), MatShape(2)))
            isScalarDia = getUnifRand()
            if (doff < 0_IK) then
                vdia = getUnifRand(vlow, vupp, min(nrow + doff, ncol))
            else
                vdia = getUnifRand(vlow, vupp, min(nrow, ncol - doff))
            end if
            if (isScalarDia) then
                call setMatInitULD_D2_ref(Mat_ref, vupp, vlow, vdia(1:), nrow, ncol, roff, coff, doff)
            else
                call setMatInitULD_D2_ref(Mat_ref, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
            end if

            call runTestsWith(doff)
            if (doff == 0_IK .and. getUnifRand()) call runTestsWith()

        end do

    contains

        subroutine runTestsWith(doff)

            integer(IK), intent(in), optional :: doff

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (allocated(mat)) deallocate(mat)
            allocate(mat(MatShape(1), MatShape(2)))
            if (roff == 0_IK .and. coff == 0_IK .and. getUnifRand()) then ! test the interfaces with NO OFFSET `roff`, `coff`.
                if (isScalarDia) then
#if                 getMatInitULD_D2_ENABLED
                    if (all([nrow, ncol] == MatShape)) then ! test for optional arguments `nrow, ncol` only when they match the matrix shape.
                        mat = getMatInit(MatShape, uppLowDia, vupp, vlow, vdia(1), doff = doff)
                    else
                        mat = getMatInit(MatShape, uppLowDia, vupp, vlow, vdia(1), nrow, ncol, doff = doff)
                    end if
#elif               setMatInitULD_D2_ENABLED
                    call setMatInit(mat, uppLowDia, vupp, vlow, vdia(1), nrow, ncol, doff, roff, coff)
#endif
                else ! vector diagonal
#if                 getMatInitULD_D2_ENABLED
                    if (all([nrow, ncol] == MatShape)) then
                        mat = getMatInit(MatShape, uppLowDia, vupp, vlow, vdia, doff = doff)
                    else
                        mat = getMatInit(MatShape, uppLowDia, vupp, vlow, vdia, nrow, ncol, doff = doff)
                    end if
#elif               setMatInitULD_D2_ENABLED
                    call setMatInit(mat, vupp, vlow, vdia, nrow, ncol, doff)
#endif
                end if
            else ! test for NON-ZERO OFFSETS for the top-left corner.
                if (isScalarDia) then ! test for scalar input diagonal.
#if                 getMatInitULD_D2_ENABLED
                    mat = getMatInit(MatShape, uppLowDia, vupp, vlow, vdia(1), nrow, ncol, roff, coff, doff)
#elif               setMatInitULD_D2_ENABLED
                    call setMatInit(mat, uppLowDia, vupp, vlow, vdia(1), nrow, ncol, roff, coff, doff)
#endif
                else ! test for vector input diagonal.
#if                 getMatInitULD_D2_ENABLED
                    mat = getMatInit(MatShape, uppLowDia, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#elif               setMatInitULD_D2_ENABLED
                    call setMatInit(mat, uppLowDia, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
#endif
                end if
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(doff)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine setMatInitULD_D2_ref(mat, vupp, vlow, vdia, nrow, ncol, roff, coff, doff)
            integer(IK), intent(in) :: nrow, ncol, roff, coff, doff
#if         getMatInitULD_D2_SK_ENABLED || setMatInitULD_D2_SK_ENABLED
            character(*,SKG) :: mat(:,:), vdia(:), vupp, vlow
#elif       getMatInitULD_D2_IK_ENABLED || setMatInitULD_D2_IK_ENABLED
            integer(IKG)     :: mat(:,:), vdia(:), vupp, vlow
#elif       getMatInitULD_D2_LK_ENABLED || setMatInitULD_D2_LK_ENABLED
            logical(LKG)     :: mat(:,:), vdia(:), vupp, vlow
#elif       getMatInitULD_D2_CK_ENABLED || setMatInitULD_D2_CK_ENABLED
            complex(CKG)     :: mat(:,:), vdia(:), vupp, vlow
#elif       getMatInitULD_D2_RK_ENABLED || setMatInitULD_D2_RK_ENABLED
            real(RKG)        :: mat(:,:), vdia(:), vupp, vlow
#else
#error      "Unrecognized Interface."
#endif
            integer(IK) :: i, irow, icol, lenDia, rowBeg, rowEnd, colBeg, colEnd
            lenDia = merge(size(vdia, 1, IK), merge(min(nrow + doff, ncol), min(nrow, ncol - doff), doff < 0), size(vdia) > 1)
            rowBeg = roff + 1_IK
            rowEnd = roff + nrow
            colBeg = coff + 1_IK
            colEnd = coff + ncol
            if (doff > 0_IK) then
                do icol = colBeg, coff + doff
                    mat(rowBeg : rowEnd, icol) = vlow
                end do
                icol = coff + doff + 1
                mat(rowBeg , icol) = vdia(1)
                mat(rowBeg + 1 : rowEnd , icol) = vlow
                do i = 2_IK, lenDia
                    irow = roff + i
                    icol = coff + i + doff
                    mat(rowBeg : irow - 1, icol) = vupp
                    mat(irow , icol) = merge(vdia(1), vdia(i), isScalarDia)
                    mat(irow + 1 : rowEnd , icol) = vlow
                end do
                do icol = coff + lenDia + doff + 1, colEnd
                    mat(rowBeg : rowEnd, icol) = vupp
                end do
            else
                do i = 1_IK, lenDia
                    irow = roff + i - doff
                    icol = coff + i
                    mat(rowBeg : irow - 1, icol) = vupp
                    mat(irow , icol) = merge(vdia(1), vdia(i), isScalarDia)
                    mat(irow + 1 : rowEnd , icol) = vlow
                end do
                do icol = coff + lenDia + 1, colEnd
                    mat(rowBeg : rowEnd, icol) = vupp
                end do
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(doff)
            use pm_io, only: display_type, field_type
            use pm_option, only: getOption
            type(display_type) :: disp
            integer(IK), optional   :: doff
            disp = display_type(test%disp%unit, deliml = field_type(string = SK_""""))
#if         getMatInitULD_D2_ENABLED
            assertion = assertion .and. logical(all(shape(mat) == shape(Mat_ref)), LK)
            call test%assert(assertion, SK_": The shape of the output matrix must conform with the reference matrix. "//getStr([shape(mat) == shape(Mat_ref)]), int(__LINE__, IK))
#endif
            assertion = assertion .and. logical(all(mat(roff + 1 : roff + nrow, coff + 1 : coff + ncol) IS_EQUAL Mat_ref(roff + 1 : roff + nrow, coff + 1 : coff + ncol)), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip()
                call disp%show("[vlow, vupp, vdia]")
                call disp%show( [vlow, vupp, vdia] )
                call disp%show("[nrow, ncol, roff, coff]")
                call disp%show( [nrow, ncol, roff, coff] )
                call disp%show("present(doff)")
                call disp%show( present(doff) )
                call disp%show("getOption(0_IK, doff)")
                call disp%show( getOption(0_IK, doff) )
                call disp%show("Mat_ref")
                call disp%show( Mat_ref )
                call disp%show("mat")
                call disp%show( mat )
                call disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_": The elements of the output matrix must match the elements of the reference matrix.", int(__LINE__, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  GET_INDEX
#undef  IS_EQUAL
#undef  GET_SIZE