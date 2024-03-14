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
!>  This include file contains procedure implementation of the generic interfaces of [pm_matrixChol](@ref pm_matrixChol).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Thursday 01:00 AM, September 23, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
        complex(TKC), parameter :: ZERO = (0._TKC, 0._TKC), ONE = (1._TKC, 0._TKC)
#define TYPE_KIND complex(TKC)
#define GET_CONJG(X) conjg(X)
#define GET_RE(x) x%re
#elif   RK_ENABLED
        real(TKC), parameter :: ZERO = 0._TKC, ONE = 1._TKC
#define TYPE_KIND real(TKC)
#define GET_CONJG(X) X
#define GET_RE(x) x
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%
#if     setChoLow_ENABLED
        !%%%%%%%%%%%%%%%%

        real(TKC) :: summ
        integer(IK) :: idim
        do idim = 1_IK, ndim
            summ = mat(idim,idim) - dot_product(mat(idim,1:idim-1), mat(idim,1:idim-1))
            if (0._TKC < summ) then
                dia(idim) = sqrt(summ)
                mat(idim+1:ndim,idim) = (mat(idim,idim+1:ndim) - matmul(mat(idim+1:ndim,1:idim-1), mat(idim,1:idim-1))) / dia(idim)
            else
                dia(1) = -idim
                return
            end if
        end do

        !%%%%%%%%%%%%%%%%%
#elif   getMatChol_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: info
#if     ULD_ENABLED
        type(uppDia_type), parameter :: subset = uppDia
#elif   !(UXD_ENABLED || XLD_ENABLED)
#error  "Unrecognized interface."
#endif
        chol = ZERO
        optionalBlock: block
            if (present(operation)) then
                if (same_type_as(operation, transHerm)) then
                    call setMatChol(mat, subset, info, chol, transHerm)
                    exit optionalBlock
                elseif (.not. same_type_as(operation, nothing)) then
                    error stop MODULE_NAME//SK_"@getMatChol(): Unrecognized `operation` other than `nothing` or `transHerm` specified." ! LCOV_EXCL_LINE
                end if
            end if
            ! default operation.
            call setMatChol(mat, subset, info, chol, nothing)
        end block optionalBlock
        if (info /= 0_IK) error stop MODULE_NAME//SK_"@getMatChol(): Cholesky factorization failed."

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatChol_ENABLED && AXX_ENABLED && ONO_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setMatChol(mat, subset, info, iteration)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatChol_ENABLED && ANI_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: summ
        integer(IK) :: irow
#if     IMP_ENABLED
        integer(IK) :: ndim
        ndim = size(mat, 1, IK)
        CHECK_ASSERTION(__LINE__, size(mat, 1, IK) == size(mat, 2, IK), SK_"@setMatChol(): The condition `size(mat, 1) == size(mat, 2)` must hold. ndim = "//getStr([size(mat, 1, IK), size(mat, 2, IK)]))
        CHECK_ASSERTION(__LINE__, all(shape(mat, IK) == shape(chol, IK)), SK_"@setMatChol(): The condition `all(shape(mat) == shape(chol))` must hold. shape(mat), shape(chol) = "//getStr([shape(mat, IK), shape(chol, IK)]))
#elif   EXP_ENABLED
        ! Ensure offsets and subset bounds are logical.
        CHECK_ASSERTION(__LINE__, 0_IK <= ndim, SK_"@setMatChol(): The condition `0 <= ndim` must hold. ndim = "//getStr(ndim))
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(mat, kind = IK)), SK_"@setMatChol(): The condition `all(0 <= [roff, coff])` must hold. roff, coff = "//getStr(1_IK - lbound(mat)))
        CHECK_ASSERTION(__LINE__, all(ndim <= ubound(mat, kind = IK)), SK_"@setMatChol(): The condition `ndim + [roff, coff] <= shape(mat)` must hold. ndim, roff, coff, shape(mat) = "//getStr([ndim, 1 - lbound(mat, kind = IK), shape(mat, kind = IK)]))
        CHECK_ASSERTION(__LINE__, all(ndim <= ubound(chol, kind = IK)), SK_"@setMatChol(): The condition `ndim + [roffc, coffc] <= shape(chol)` must hold. ndim, roffc, coffc, shape(chol) = "//getStr([ndim, 1 - lbound(chol, kind = IK), shape(chol, kind = IK)]))
#else
#error  "Unrecognized interface."
#endif
        ! Define operations based on the specified triangular side.
#if     UXD_ENABLED
#define GET_MATMUL(I,J) matmul(J,I)
#define GETMAT(MAT,I,J) MAT(J,I)
#elif   XLD_ENABLED
#define GET_MATMUL(I,J) matmul(I,J)
#define GETMAT(MAT,I,J) MAT(I,J)
#else
#error  "Unrecognized interface."
#endif
        do info = 1_IK, ndim
#if         ONO_ENABLED
            summ = real(mat(info, info), TKC) - real(dot_product(GETMAT(chol, info, 1 : info - 1), GETMAT(chol, info, 1 : info - 1)), TKC)
#elif       OTH_ENABLED
            summ = real(mat(info, info), TKC) - real(dot_product(GETMAT(chol, 1 : info - 1, info), GETMAT(chol, 1 : info - 1, info)), TKC)
#else
#error      "Unrecognized interface."
#endif
            if (0._TKC < summ) then
                summ = sqrt(summ)
                chol(info, info) = summ
                summ = 1._TKC / summ
#if             ONO_ENABLED
                !GETMAT(chol, info + 1 : ndim, info) = (GETMAT(mat, info + 1 : ndim, info) - GET_MATMUL(GET_CONJG(GETMAT(chol, info + 1 : ndim, 1 : info - 1)),GETMAT(chol, info, 1 : info - 1))) * summ
                do irow = info + 1, ndim
                    GETMAT(chol, irow, info) = (GETMAT(mat, irow, info) - dot_product(GETMAT(chol, info, 1 : info - 1), GETMAT(chol, irow, 1 : info - 1))) * summ
                end do
#elif           OTH_ENABLED
                do irow = info + 1, ndim
                    GETMAT(chol, info, irow) = (GET_CONJG(GETMAT(mat, irow, info) - dot_product(GETMAT(chol, 1 : info - 1, irow), GETMAT(chol, 1 : info - 1, info)))) * summ
                end do
#endif
                cycle
            end if
            return
        end do
        info = 0_IK
#undef  GET_MATMUL
#undef  GETMAT

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatChol_ENABLED && (ABI_ENABLED || ABR_ENABLED) && ONO_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define OPERATION transHerm
#elif   RK_ENABLED
#define OPERATION transSymm
#else
#error  "Unrecognized interface."
#endif
#if     ABI_ENABLED
        ! Set the block size for this environment.
        ! Per LAPACK documentation, the current value (ILAENV( 1, 'dpotrf', uplo, ndim, -1, -1, -1 )) may not be optimal.
        ! Further runtime experimentation might be needed.
        integer(IK) , parameter :: NBLOCK = 64_IK
        integer(IK) :: icol, jcol, bdim_def
#elif   ABR_ENABLED
        integer(IK) :: ndimHalf1, ndimHalf2
#endif
        ! Set/check bounds.
#if     IMP_ENABLED
        integer(IK) , parameter :: roff = 0, coff = 0
        integer(IK) :: ndim
        ndim = size(mat, 1, IK)
        CHECK_ASSERTION(__LINE__, size(mat, 1, IK) == size(mat, 2, IK), SK_"@setMatChol(): The condition `size(mat, 1) == size(mat, 2)` must hold. ndim = "//getStr([size(mat, 1, IK), size(mat, 2, IK)]))
#elif   EXP_ENABLED
        ! Ensure offsets and subset bounds are logical.
        CHECK_ASSERTION(__LINE__, 0_IK <= ndim, SK_"@setMatChol(): The condition `0 <= ndim` must hold. ndim = "//getStr(ndim))
        CHECK_ASSERTION(__LINE__, all(0_IK <= 1_IK - lbound(mat, kind = IK)), SK_"@setMatChol(): The condition `all(0 <= [roff, coff])` must hold. roff, coff = "//getStr(1_IK - lbound(mat)))
        CHECK_ASSERTION(__LINE__, all(ndim <= ubound(mat, kind = IK)), SK_"@setMatChol(): The condition `ndim + [roff, coff] <= shape(mat)` must hold. ndim, roff, coff, shape(mat) = "//getStr([ndim, 1 - lbound(mat, kind = IK), shape(mat, kind = IK)]))
#else
#error  "Unrecognized interface."
#endif
        info = 0_IK
        if (ndim == 0_IK) return ! Quick return if possible.

        !%%%%%%%%%%
#if     ABI_ENABLED
        !%%%%%%%%%%

        if (present(bdim)) then
            CHECK_ASSERTION(__LINE__, 0_IK < bdim, SK_"@setMatChol(): The condition `0 < bdim` must hold. bdim = "//getStr(bdim))
            bdim_def = bdim
        else
            bdim_def = NBLOCK
        end if
        if (bdim_def <= 1_IK .or. ndim <= bdim_def) then ! Use unblocked code.
            call setMatChol(mat, subset, info, recursion, ndim, roff, coff)
        else ! Use blocked code.
            ! Compute the Cholesky factorization mat = U**T*U.
            do icol = 0_IK, ndim - 1_IK, bdim_def
                ! Update and factorize the current diagonal block and test for non-positive-definiteness.
                jcol = min(bdim_def, ndim - icol)
#if             UXD_ENABLED
                call setMatUpdate   ( mat = mat & ! LCOV_EXCL_LINE
                                    , class = hermitian & ! LCOV_EXCL_LINE
                                    , subset = uppDia & ! LCOV_EXCL_LINE
                                    , matA = mat & ! LCOV_EXCL_LINE
                                    , operationA = trans & ! LCOV_EXCL_LINE
                                    , alpha = -ONE & ! LCOV_EXCL_LINE
                                    , beta = ONE & ! LCOV_EXCL_LINE
                                    , ndim = jcol & ! LCOV_EXCL_LINE
                                    , ndum = icol & ! LCOV_EXCL_LINE
                                    , roff = roff + icol & ! LCOV_EXCL_LINE
                                    , coff = coff + icol & ! LCOV_EXCL_LINE
                                    , roffA = roff & ! LCOV_EXCL_LINE
                                    , coffA = coff + icol & ! LCOV_EXCL_LINE
                                    ) ! Update and factor MatPosDef22.
                call setMatChol(mat, subset, info, recursion, jcol, roff + icol, coff + icol)
                if (info /= 0_IK) then
                    info = info + icol
                    return
                end if
                if (icol + jcol <= ndim) then ! Compute the current block row.
                    call setMatMulAdd   ( matA = mat & ! LCOV_EXCL_LINE
                                        , operationA = transHerm & ! LCOV_EXCL_LINE
                                        , matB = mat & ! LCOV_EXCL_LINE
                                        , matC = mat & ! LCOV_EXCL_LINE
                                        , alpha = -ONE & ! LCOV_EXCL_LINE
                                        , beta = ONE & ! LCOV_EXCL_LINE
                                        , nrow = jcol & ! LCOV_EXCL_LINE
                                        , ncol = ndim - icol - jcol & ! LCOV_EXCL_LINE
                                        , ndum = icol & ! LCOV_EXCL_LINE
                                        , roffA = roff & ! LCOV_EXCL_LINE
                                        , coffA = coff + icol & ! LCOV_EXCL_LINE
                                        , roffB = roff & ! LCOV_EXCL_LINE
                                        , coffB = coff + icol + jcol & ! LCOV_EXCL_LINE
                                        , roffC = roff + icol & ! LCOV_EXCL_LINE
                                        , coffC = coff + icol + jcol & ! LCOV_EXCL_LINE
                                        )
                    ! syslin = AXB
                    call setMatMulTri   ( matA = mat & ! LCOV_EXCL_LINE
                                        , classA = upperDiag & ! LCOV_EXCL_LINE
                                        , operationA = transUnit & ! LCOV_EXCL_LINE
                                        , matB = mat & ! LCOV_EXCL_LINE
                                        , alpha = ONE & ! LCOV_EXCL_LINE
                                        , nrow = jcol & ! LCOV_EXCL_LINE
                                        , ncol = ndim - icol - jcol & ! LCOV_EXCL_LINE
                                        , roffA = roff + icol & ! LCOV_EXCL_LINE
                                        , coffA = coff + icol & ! LCOV_EXCL_LINE
                                        , roffB = roff + icol & ! LCOV_EXCL_LINE
                                        , coffB = coff + icol + jcol & ! LCOV_EXCL_LINE
                                        )
                end if
#elif           XLD_ENABLED
                ! Compute the Cholesky factorization mat = l*l**t.
                call setMatUpdate   ( mat = mat & ! LCOV_EXCL_LINE
                                    , class = hermitian & ! LCOV_EXCL_LINE
                                    , subset = lowDia & ! LCOV_EXCL_LINE
                                    , matA = mat & ! LCOV_EXCL_LINE
                                    , operationA = nothing & ! LCOV_EXCL_LINE
                                    , alpha = -ONE & ! LCOV_EXCL_LINE
                                    , beta = ONE & ! LCOV_EXCL_LINE
                                    , ndim = jcol & ! LCOV_EXCL_LINE
                                    , ndum = icol & ! LCOV_EXCL_LINE
                                    , roff = roff + icol & ! LCOV_EXCL_LINE
                                    , coff = coff + icol & ! LCOV_EXCL_LINE
                                    , roffA = roff + icol & ! LCOV_EXCL_LINE
                                    , coffA = coff & ! LCOV_EXCL_LINE
                                    ) ! Update and factor MatPosDef22.
                call setMatChol(mat, subset, info, recursion, jcol, roff + icol, coff + icol)
                if (info /= 0_IK) then
                    info = info + icol
                    return
                end if
                if (icol + jcol <= ndim) then ! Compute the current block column.
                    call setMatMulAdd   ( matA = mat & ! LCOV_EXCL_LINE
                                        , matB = mat & ! LCOV_EXCL_LINE
                                        , matC = mat & ! LCOV_EXCL_LINE
                                        , alpha = -ONE & ! LCOV_EXCL_LINE
                                        , beta = ONE & ! LCOV_EXCL_LINE
                                        , nrow = ndim - icol - jcol & ! LCOV_EXCL_LINE
                                        , ncol = jcol & ! LCOV_EXCL_LINE
                                        , ndum = icol & ! LCOV_EXCL_LINE
                                        , roffA = roff + icol + jcol & ! LCOV_EXCL_LINE
                                        , coffA = coff & ! LCOV_EXCL_LINE
                                        , roffB = roff + icol & ! LCOV_EXCL_LINE
                                        , coffB = coff & ! LCOV_EXCL_LINE
                                        , roffC = roff + icol + jcol & ! LCOV_EXCL_LINE
                                        , coffC = coff + icol & ! LCOV_EXCL_LINE
                                        , operationB = transHerm & ! LCOV_EXCL_LINE
                                        )
                    ! syslin = XAB
                    call setMatMulTri   ( matA = mat & ! LCOV_EXCL_LINE
                                        , matB = mat & ! LCOV_EXCL_LINE
                                        , classB = lowerDiag & ! LCOV_EXCL_LINE
                                        , operationB = transUnit & ! LCOV_EXCL_LINE
                                        , alpha = ONE & ! LCOV_EXCL_LINE
                                        , nrow = ndim - icol - jcol & ! LCOV_EXCL_LINE
                                        , ncol = jcol & ! LCOV_EXCL_LINE
                                        , roffA = roff + icol + jcol & ! LCOV_EXCL_LINE
                                        , coffA = coff + icol & ! LCOV_EXCL_LINE
                                        , roffB = roff + icol & ! LCOV_EXCL_LINE
                                        , coffB = coff + icol & ! LCOV_EXCL_LINE
                                        )
                end if
#else
#error          "Unrecognized interface."
#endif
            end do
        end if

        !%%%%%%%%%%
#elif   ABR_ENABLED
        !%%%%%%%%%%

        if(1_IK < ndim) then ! use recursive code
            ndimHalf1 = ndim / 2_IK
            ndimHalf2 = ndim - ndimHalf1
            call setMatChol(mat, subset, info, control, ndimHalf1, roff, coff) ! Factor MatPosDef11.
            if (info /= 0_IK) return
#if         UXD_ENABLED
            ! Compute the Cholesky factorization mat = U**T*U :: syslin = AXB
            call setMatMulTri   ( matA = mat & ! LCOV_EXCL_LINE
                                , classA = upperDiag & ! LCOV_EXCL_LINE
                                , operationA = transUnit & ! LCOV_EXCL_LINE
                                , matB = mat & ! LCOV_EXCL_LINE
                                , alpha = ONE & ! LCOV_EXCL_LINE
                                , nrow = ndimHalf1 & ! LCOV_EXCL_LINE
                                , ncol = ndimHalf2 & ! LCOV_EXCL_LINE
                                , roffA = roff & ! LCOV_EXCL_LINE
                                , coffA = coff & ! LCOV_EXCL_LINE
                                , roffB = roff & ! LCOV_EXCL_LINE
                                , coffB = coff + ndimHalf1 & ! LCOV_EXCL_LINE
                                ) ! Update and scale MatPosDef12.
            !block
            !use pm_blas, only: blasTRSM
            !CALL blasTRSM('L', 'U', 'C', 'N', ndimHalf1, ndimHalf2, ONE, mat(1, 1), size(mat, 1), mat(1, ndimHalf1 + 1), size(mat, 1))
            !end block
            call setMatUpdate   ( mat = mat & ! LCOV_EXCL_LINE
                                , class = hermitian & ! LCOV_EXCL_LINE
                                , subset = uppDia & ! LCOV_EXCL_LINE
                                , matA = mat & ! LCOV_EXCL_LINE
                                , operationA = trans & ! LCOV_EXCL_LINE
                                , alpha = -ONE & ! LCOV_EXCL_LINE
                                , beta = ONE & ! LCOV_EXCL_LINE
                                , ndim = ndimHalf2 & ! LCOV_EXCL_LINE
                                , ndum = ndimHalf1 & ! LCOV_EXCL_LINE
                                , roff = roff + ndimHalf1 & ! LCOV_EXCL_LINE
                                , coff = coff + ndimHalf1 & ! LCOV_EXCL_LINE
                                , roffA = roff & ! LCOV_EXCL_LINE
                                , coffA = coff + ndimHalf1 & ! LCOV_EXCL_LINE
                                ) ! Update and factor MatPosDef22.
#elif       XLD_ENABLED
            ! Compute the Cholesky factorization mat = L*L**T :: syslin = XAB
            call setMatMulTri   ( matA = mat & ! LCOV_EXCL_LINE
                                , matB = mat & ! LCOV_EXCL_LINE
                                , classB = lowerDiag & ! LCOV_EXCL_LINE
                                , operationB = transUnit & ! LCOV_EXCL_LINE
                                , alpha = ONE & ! LCOV_EXCL_LINE
                                , nrow = ndimHalf2 & ! LCOV_EXCL_LINE
                                , ncol = ndimHalf1 & ! LCOV_EXCL_LINE
                                , roffA = roff + ndimHalf1 & ! LCOV_EXCL_LINE
                                , coffA = coff & ! LCOV_EXCL_LINE
                                , roffB = roff & ! LCOV_EXCL_LINE
                                , coffB = coff & ! LCOV_EXCL_LINE
                                ) ! Update and scale MatPosDef21.
            call setMatUpdate   ( mat = mat & ! LCOV_EXCL_LINE
                                , class = hermitian & ! LCOV_EXCL_LINE
                                , subset = lowDia & ! LCOV_EXCL_LINE
                                , matA = mat & ! LCOV_EXCL_LINE
                                , operationA = nothing & ! LCOV_EXCL_LINE
                                , alpha = -ONE & ! LCOV_EXCL_LINE
                                , beta = ONE & ! LCOV_EXCL_LINE
                                , ndim = ndimHalf2 & ! LCOV_EXCL_LINE
                                , ndum = ndimHalf1 & ! LCOV_EXCL_LINE
                                , roff = roff + ndimHalf1 & ! LCOV_EXCL_LINE
                                , coff = coff + ndimHalf1 & ! LCOV_EXCL_LINE
                                , roffA = roff + ndimHalf1 & ! LCOV_EXCL_LINE
                                , coffA = coff & ! LCOV_EXCL_LINE
                                ) ! Update and factor MatPosDef22.
#else
#error      "Unrecognized interface."
#endif
            call setMatChol(mat, subset, info, control, ndimHalf2, roff + ndimHalf1, coff + ndimHalf1)
            if (info /= 0_IK) then
                info = info + ndimHalf1
                return
            end if
        else
            ! Test for non-positive-definiteness.
            if (GET_RE(ZERO) < GET_RE(mat(1, 1))) then
                mat(1, 1) = sqrt(GET_RE(mat(1, 1)))
            else
                !write(*,*) GET_RE(mat(1, 1))
                info = 1_IK
                return
            end if
        end if
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_KIND
#undef  GET_CONJG
#undef  OPERATION
#undef  GETMAT
#undef  GET_RE
#undef  CHOMAT