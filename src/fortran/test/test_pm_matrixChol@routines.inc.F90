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
!>  This include file contains the implementations of the tests of procedures of [pm_matrixChol](@ref pm_matrixChol).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG), parameter :: RTOL = sqrt(epsilon(0._TKG))
#if     CK_ENABLED
#define TYPE_OF cmplx
#define GET_CONJG(X) conjg(X)
#define TYPE_OF_MAT complex(TKG)
        complex(TKG), parameter :: ZERO = (0._TKG, 0._TKG), ONE = (1._TKG, 0._TKG), LB = (1._TKG, -1._TKG), UB = (2._TKG, 1._TKG), TOL = (RTOL, RTOL)
#elif   RK_ENABLED
#define TYPE_OF real
#define GET_CONJG(X) X
#define TYPE_OF_MAT real(TKG)
        real(TKG), parameter :: ZERO = 0._TKG, ONE = 1._TKG, LB = 1._TKG, UB = 2._TKG, TOL = RTOL
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%
#if     setChoLow_ENABLED
        !%%%%%%%%%%%%%%%%

        integer(IK), parameter :: ntry = 100_IK
        integer(IK) :: ndim, info, info_ref, itry
        TYPE_OF_MAT, allocatable :: chol_ref(:,:), chol(:,:), mat(:,:), diff(:,:), vdia(:), vdia_ref(:)

        assertion = .true._LK

        do itry = 1, ntry

            ! Set the matrix rank.
            ndim = getUnifRand(1_IK, 7_IK)
            chol_ref = getFilled(ZERO, ndim, ndim)
            ! Generate random upper-triangular matrix.
            mat = getCovRand(chol_ref(1,1), ndim)
            call setMatInit(mat, low, ZERO)
            vdia = getFilled(ZERO, ndim)
            info = 0_IK

            call setMatChol(mat, uppDia, info_ref, chol_ref, transHerm)
            if (info_ref /= 0_IK) error stop getFine(__FILE__, __LINE__)//SK_": setMatChol() failed." ! LCOV_EXCL_LINE
            vdia_ref = getMatCopy(lfpack, chol_ref, rdpack, dia)
            call setMatCopy(chol_ref, rdpack, mat, rdpack, uppDia)
            chol = mat
            call setChoLow(chol, vdia, ndim)
            diff = chol - chol_ref
            assertion = assertion .and. all(diff < TOL)
            assertion = assertion .and. all(abs(vdia - vdia_ref) < TOL)
            call report(__LINE__)

            info_ref = getUnifRand(1_IK, ndim)
            mat(info_ref, info_ref) = -mat(info_ref, info_ref)
            call setMatChol(mat, uppDia, info_ref, chol_ref, transHerm)
            vdia_ref = getMatCopy(lfpack, chol_ref, rdpack, dia)
            call setMatCopy(chol_ref, rdpack, mat, rdpack, uppDia)
            chol = mat
            call setChoLow(chol, vdia, ndim)
            info = int(-vdia(1), IK)
            assertion = assertion .and. info == info_ref
            call report(__LINE__)

        end do

    contains

        subroutine report(line)
            integer, intent(in) :: line
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip
                call test%disp%show("ndim")
                call test%disp%show( ndim )
                call test%disp%show("mat")
                call test%disp%show( mat )
                call test%disp%show("chol_ref")
                call test%disp%show( chol_ref )
                call test%disp%show("chol")
                call test%disp%show( chol )
                call test%disp%show("diff")
                call test%disp%show( diff )
                call test%disp%show("vdia_ref")
                call test%disp%show( vdia_ref )
                call test%disp%show("vdia")
                call test%disp%show( vdia )
                call test%disp%show("vdia - vdia_ref")
                call test%disp%show( vdia - vdia_ref )
                call test%disp%show("info_ref")
                call test%disp%show( info_ref )
                call test%disp%show("info")
                call test%disp%show( info )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The Cholesky factorization must not fail.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%
#elif   getMatChol_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim, info, itry
        integer(IK), parameter :: ntry = 100_IK
        TYPE_OF_MAT, allocatable :: chol_ref(:,:), chol(:,:), mat(:,:), diff(:,:)

        assertion = .true._LK

        do itry = 1, ntry

            ! Set the matrix rank.
            ndim = getUnifRand(1_IK, 7_IK)
            chol_ref = getFilled(ZERO, ndim, ndim)
            ! Generate random lower-triangular matrix.
            mat = getCovRand(chol_ref(1,1), ndim)
            call setMatInit(mat, upp, ZERO)

            ! lowDia

            block

                chol_ref = ZERO
                call setMatChol(mat, lowDia, info, chol_ref, nothing)
                if (info /= 0_IK) error stop getFine(__FILE__, __LINE__)//SK_": setMatChol() failed." ! LCOV_EXCL_LINE
                chol = getMatChol(mat, lowDia)
                call report(__LINE__, "lowDia")

                chol_ref = ZERO
                call setMatChol(mat, lowDia, info, chol_ref, nothing)
                if (info /= 0_IK) error stop getFine(__FILE__, __LINE__)//SK_": setMatChol() failed." ! LCOV_EXCL_LINE
                chol = getMatChol(mat, lowDia, nothing)
                call report(__LINE__, "lowDia", "nothing")

                chol_ref = ZERO
                call setMatChol(mat, lowDia, info, chol_ref, transHerm)
                if (info /= 0_IK) error stop getFine(__FILE__, __LINE__)//SK_": setMatChol() failed." ! LCOV_EXCL_LINE
                chol = getMatChol(mat, lowDia, transHerm)
                call report(__LINE__, "lowDia", "transHerm")

            end block

            mat = transpose(GET_CONJG(mat))

            ! uppDia

            block

                chol_ref = ZERO
                call setMatChol(mat, uppDia, info, chol_ref, nothing)
                if (info /= 0_IK) error stop getFine(__FILE__, __LINE__)//SK_": setMatChol() failed." ! LCOV_EXCL_LINE
                chol = getMatChol(mat, uppDia)
                call report(__LINE__, "uppDia")

                chol_ref = ZERO
                call setMatChol(mat, uppDia, info, chol_ref, nothing)
                if (info /= 0_IK) error stop getFine(__FILE__, __LINE__)//SK_": setMatChol() failed." ! LCOV_EXCL_LINE
                chol = getMatChol(mat, uppDia, nothing)
                call report(__LINE__, "uppDia", "nothing")

                chol_ref = ZERO
                call setMatChol(mat, uppDia, info, chol_ref, transHerm)
                if (info /= 0_IK) error stop getFine(__FILE__, __LINE__)//SK_": setMatChol() failed." ! LCOV_EXCL_LINE
                chol = getMatChol(mat, uppDia, transHerm)
                call report(__LINE__, "uppDia", "transHerm")

            end block

        end do

    contains

        subroutine report(line, subset, operation)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: subset
            character(*, SK), intent(in), optional :: operation
            diff = chol - chol_ref
            assertion = assertion .and. all(diff < TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip
                call test%disp%show("ndim")
                call test%disp%show( ndim )
                call test%disp%show("mat")
                call test%disp%show( mat )
                call test%disp%show("subset")
                call test%disp%show( subset )
                call test%disp%show("present(operation)")
                call test%disp%show( present(operation) )
                if (present(operation)) then
                    call test%disp%show("operation")
                    call test%disp%show( operation )
                end if
                call test%disp%show("chol_ref")
                call test%disp%show( chol_ref )
                call test%disp%show("chol")
                call test%disp%show( chol )
                call test%disp%show("diff")
                call test%disp%show( diff )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The Cholesky factorization must not fail.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%
#elif   setMatChol_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim, info, info_def, itry
        integer(IK), parameter :: ntry = 100_IK
        TYPE_OF_MAT, allocatable :: choUpp_ref(:,:), choLow_ref(:,:), mat_ref(:,:)
        TYPE_OF_MAT, allocatable :: matUpp(:,:), matLow(:,:), choUpp(:,:), choLow(:,:), diff(:,:)

        info_def = 0_IK
        assertion = .true._LK

        do itry = 1, ntry

            ! Set the matrix rank.
            ndim = getUnifRand(1_IK, 7_IK)
            ! Generate random triangular matrix.
            choLow_ref = getUnifRand(LB, UB, ndim, ndim)
            choLow_ref = getMatCopy(rdpack, choLow_ref, rdpack, lowDia, init = ZERO)
            ! set the diagonals to positive real values.
            call setMatInit(choLow_ref, dia, TYPE_OF(getUnifRand(1._TKG, 2._TKG, ndim), kind = TKG))
            ! Create the reference upper Cholesky matrix.
            choUpp_ref = transpose(GET_CONJG(choLow_ref))
            ! Create the positive definite matrix.
            mat_ref = matmul(choLow_ref, choUpp_ref)
            matUpp = getMatCopy(rdpack, mat_ref, rdpack, uppDia, init = ZERO)
            matLow = getMatCopy(rdpack, mat_ref, rdpack, lowDia, init = ZERO)

            ! Test the routines.

            internal_interface: block

                ! internal paramonte interface, no overwrite.

                choUpp = getFilled(ZERO, ndim, ndim)
                choLow = getFilled(ZERO, ndim, ndim)

                choUpp = getFilled(ZERO, ndim, ndim)
                call setMatChol(matUpp, uppDia, info, choUpp, nothing)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "uppDia", "nothing", .false._LK, info)
                diff = abs(choUpp - choUpp_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "uppDia", "nothing", .false._LK, diff = diff)

                choLow = getFilled(ZERO, ndim, ndim)
                call setMatChol(matLow, lowDia, info, choLow, nothing)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "lowDia", "nothing", .false._LK, info)
                diff = abs(choLow - choLow_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "lowDia", "nothing", .false._LK, diff = diff)

                choLow = getFilled(ZERO, ndim, ndim)
                call setMatChol(matUpp, uppDia, info, choLow, transHerm)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "uppDia", "tranHerm", .false._LK, info)
                diff = abs(choLow - choLow_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "uppDia", "tranHerm", .false._LK, diff = diff)

                choUpp = getFilled(ZERO, ndim, ndim)
                call setMatChol(matLow, lowDia, info, choUpp, transHerm)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "lowDia", "transHerm", .false._LK, info)
                diff = abs(choUpp - choUpp_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "lowDia", "transHerm", .false._LK, diff = diff)

                ! internal paramonte interface, with overwrite.

                choUpp = getFilled(ZERO, ndim, ndim + 1)
                choLow = getFilled(ZERO, ndim, ndim + 1)

                choUpp(:,1:ndim) = matUpp
                call setMatChol(choUpp(:,1:ndim), uppDia, info, choUpp(:,1:ndim), nothing)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "uppDia", "nothing", .true._LK, info)
                diff = getDiff(choUpp(:,1:ndim), choUpp_ref, uppDia)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "uppDia", "nothing", .true._LK, diff = diff)

                choLow(:,1:ndim) = matLow
                call setMatChol(choLow(:,1:ndim), lowDia, info, choLow(:,1:ndim), nothing)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "lowDia", "nothing", .true._LK, info)
                diff = getDiff(choLow(:,1:ndim), choLow_ref, lowDia)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "lowDia", "nothing", .true._LK, diff = diff)

                choUpp(:,1:ndim) = matLow
                call setMatChol(choUpp(:,1:ndim), lowDia, info, choUpp(:,2:ndim+1), transHerm)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "uppDia", "transHerm", .true._LK, info)
                diff = getDiff(choUpp(:,2:ndim+1), choUpp_ref, uppDia)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "uppDia", "transHerm", .true._LK, diff = diff)

                choLow(:,2:ndim+1) = matUpp
                call setMatChol(choLow(:,2:ndim+1), uppDia, info, choLow(:,1:ndim), transHerm)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "lowDia", "transHerm", .true._LK, info)
                diff = getDiff(choLow(:,1:ndim), choLow_ref, lowDia)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "lowDia", "transHerm", .true._LK, diff = diff)

            end block internal_interface

            recursion_interface: block

                ! implicit interface.

                choUpp = getFilled(ZERO, ndim, ndim)
                choLow = getFilled(ZERO, ndim, ndim)

                choUpp = getMatCopy(rdpack, matUpp, rdpack, uppDia, init = ZERO)
                call setMatChol(choUpp, uppDia, info, recursion)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "uppDia", "recursion", .true._LK, info)
                diff = abs(choUpp - choUpp_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "uppDia", "recursion", .true._LK, diff = diff)

                choLow = getMatCopy(rdpack, matLow, rdpack, lowDia, init = ZERO)
                call setMatChol(choLow, lowDia, info, recursion)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "lowDia", "recursion", .true._LK, info)
                diff = abs(choLow - choLow_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "lowDia", "recursion", .true._LK, diff = diff)

            end block recursion_interface

            iteration_interface: block

                ! implicit interface.

                choUpp = getFilled(ZERO, ndim, ndim)
                choLow = getFilled(ZERO, ndim, ndim)

                choUpp = getMatCopy(rdpack, matUpp, rdpack, uppDia, init = ZERO)
                call setMatChol(choUpp, uppDia, info, iteration)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "uppDia", "iteration", .true._LK, info)
                diff = abs(choUpp - choUpp_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "uppDia", "iteration", .true._LK, diff = diff)

                choLow = getMatCopy(rdpack, matLow, rdpack, lowDia, init = ZERO)
                call setMatChol(choLow, lowDia, info, iteration)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "lowDia", "iteration", .true._LK, info)
                diff = abs(choLow - choLow_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "lowDia", "iteration", .true._LK, diff = diff)

            end block iteration_interface

            iteration_bdim_interface: block

                ! implicit interface.

                integer(IK) :: bdim

                bdim = getUnifRand(2_IK, 2_IK * ndim + 1_IK)
                choUpp = getFilled(ZERO, ndim, ndim)
                choLow = getFilled(ZERO, ndim, ndim)

                choUpp = getMatCopy(rdpack, matUpp, rdpack, uppDia, init = ZERO)
                call setMatChol(choUpp, uppDia, info, iteration, bdim = bdim)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "uppDia", "iteration", .true._LK, info, bdim = bdim)
                diff = abs(choUpp - choUpp_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "uppDia", "iteration", .true._LK, diff = diff, bdim = bdim)

                choLow = getMatCopy(rdpack, matLow, rdpack, lowDia, init = ZERO)
                call setMatChol(choLow, lowDia, info, iteration, bdim = bdim)
                assertion = assertion .and. info == 0_IK
                call report(__LINE__, "lowDia", "iteration", .true._LK, info, bdim = bdim)
                diff = abs(choLow - choLow_ref)
                assertion = assertion .and. all(diff < TOL)
                call report(__LINE__, "lowDia", "iteration", .true._LK, diff = diff, bdim = bdim)

            end block iteration_bdim_interface

            nonposdef_interface: block

                ! implicit interface.

                integer(IK) :: bdim

                info_def = getUnifRand(1_IK, ndim)
                bdim = getUnifRand(2_IK, 2_IK * ndim + 1_IK)
                mat_ref = getMatInit([ndim, ndim], uppLowDia, ZERO, ZERO, getUnifRand(LB, UB, ndim))
                mat_ref(info_def, info_def) = -ONE

                ! default

                choUpp = mat_ref
                call setMatChol(mat_ref, uppDia, info, choUpp, nothing)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "uppDia", "nothing", .false._LK, info)

                choLow = mat_ref
                call setMatChol(mat_ref, lowDia, info, choLow, nothing)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "lowDia", "nothing", .false._LK, info)

                choUpp = mat_ref
                call setMatChol(mat_ref, uppDia, info, choUpp, transHerm)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "uppDia", "transHerm", .false._LK, info)

                choLow = mat_ref
                call setMatChol(mat_ref, lowDia, info, choLow, transHerm)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "lowDia", "transHerm", .false._LK, info)

                choUpp = mat_ref
                call setMatChol(choUpp, uppDia, info)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "uppDia", "nothing", .true._LK, info)

                choLow = mat_ref
                call setMatChol(choLow, lowDia, info)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "lowDia", "nothing", .true._LK, info)

                ! recursion

                choUpp = mat_ref
                call setMatChol(choUpp, uppDia, info, recursion)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "uppDia", "recursion", .true._LK, info)

                choLow = mat_ref
                call setMatChol(choLow, lowDia, info, recursion)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "lowDia", "recursion", .true._LK, info)

                ! iteration

                choUpp = mat_ref
                call setMatChol(choUpp, uppDia, info, iteration)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "uppDia", "iteration", .true._LK, info)

                choLow = mat_ref
                call setMatChol(choLow, lowDia, info, iteration)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "lowDia", "iteration", .true._LK, info)

                choUpp = mat_ref
                call setMatChol(choUpp, uppDia, info, iteration, bdim = bdim)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "uppDia", "iteration", .true._LK, info, bdim = bdim)

                choLow = mat_ref
                call setMatChol(choLow, lowDia, info, iteration, bdim = bdim)
                assertion = assertion .and. info == info_def
                call report(__LINE__, "lowDia", "iteration", .true._LK, info, bdim = bdim)

            end block nonposdef_interface

        end do

    contains

        pure function getDiff(mat, ref, subset) result(diff)
            TYPE_OF_MAT, intent(in) :: mat(:,:), ref(:,:)
            TYPE_OF_MAT :: diff(size(ref, 1, IK), size(ref, 2, IK))
            class(subset_type), intent(in) :: subset
            integer(IK) :: idim
            diff = 0._TKG
            do idim = 1, size(ref, 1, IK)
                if (same_type_as(subset, uppDia)) then
                    diff(idim, idim :) = abs(mat(idim, idim :) - ref(idim, idim :))
                elseif (same_type_as(subset, lowDia)) then
                    diff(idim :, idim) = abs(mat(idim :, idim) - ref(idim :, idim))
                else
                    error stop "Internal library error: Unrecognized `subset` type. Please report this error to the developers."
                end if
            end do
        end function

        subroutine report(line, subset, operation, overwrite, info, diff, bdim)
            integer, intent(in) :: line
            integer, intent(in), optional :: info
            integer, intent(in), optional :: bdim
            character(*, SK), intent(in) :: subset, operation
            TYPE_OF_MAT, intent(in), optional :: diff(:,:)
            logical(LK), intent(in) :: overwrite
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip
                call test%disp%show("ndim")
                call test%disp%show( ndim )
                call test%disp%show("choUpp_ref")
                call test%disp%show( choUpp_ref )
                call test%disp%show("choUpp")
                call test%disp%show( choUpp )
                call test%disp%show("choUpp - choUpp_ref")
                call test%disp%show( choUpp - choUpp_ref )
                call test%disp%show("choLow_ref")
                call test%disp%show( choLow_ref )
                call test%disp%show("choLow")
                call test%disp%show( choLow )
                call test%disp%show("choLow - choLow_ref")
                call test%disp%show( choLow - choLow_ref )
                call test%disp%show("mat_ref")
                call test%disp%show( mat_ref )
                call test%disp%show("subset")
                call test%disp%show( subset )
                call test%disp%show("operation")
                call test%disp%show( operation )
                call test%disp%show("overwrite")
                call test%disp%show( overwrite )
                call test%disp%show("present(info)")
                call test%disp%show( present(info) )
                if (present(info)) then
                    call test%disp%show("info")
                    call test%disp%show( info )
                end if
                call test%disp%show("info_def")
                call test%disp%show( info_def )
                if (present(diff)) then
                    call test%disp%show("diff")
                    call test%disp%show( diff )
                end if
                call test%disp%show("present(bdim)")
                call test%disp%show( present(bdim) )
                if (present(bdim)) then
                    call test%disp%show("bdim")
                    call test%disp%show( bdim )
                end if
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The Cholesky factorization must not fail.", int(line, IK))
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interrface."
        !%%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_MAT
#undef  GET_CONJG
#undef  TYPE_OF