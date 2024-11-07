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
!>  This include file contains procedure implementation of the generic interface [pm_matrixLUP](@ref pm_matrixLUP).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !real(TKG), parameter :: IHUGE = 1._TKG / huge(0._TKG)
        !real(TKG), parameter :: SMALL = merge(tiny(0._TKG), IHUGE * (1._TKG + epsilon(0._TKG)), IHUGE < tiny(0._TKG))
        real(TKG), parameter :: SMALL = tiny(0._TKG)**.999
#if     CK_ENABLED
#define GET_CONJG(X) conjg(X)
#define TYPE_KIND complex(TKG)
#define GET_ABSL1(X) abs(X%re) + abs(X%im)
        complex(TKG), parameter :: ZERO = (0._TKG, 0._TKG), ONE = (1._TKG, 0._TKG)
#elif   RK_ENABLED
#define GET_CONJG(X) X
#define GET_ABSL1(X) abs(X)
#define TYPE_KIND real(TKG)
        real(TKG), parameter :: ZERO = 0._TKG, ONE = 1._TKG
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     setMatLUP_ENABLED && SQM_ENABLED && LAPACK_ENABLED && DISPATCH_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !check_assertion(__LINE__, all(size(rperm, 1, IK) == shape(mat, IK)), \
        !SK_"@setMatLUP(): The condition `all(size(rperm) == shape(mat))` must hold. size(rperm), shape(mat) = "//\
        !getStr([size(rperm, 1, IK), shape(mat, IK)])) ! fpp
        call lapackGETRF(size(mat, 1, IK), size(mat, 2, IK), mat(1,1), size(mat, 1, IK), rperm, info)
        !if (present(parity)) then
        !    block
        !        integer(IK) :: idim
        !        parity = ONE
        !        do idim = 1, size(rperm, 1, IK)
        !            if (rperm(idim) /= idim) parity = -parity
        !        end do
        !    end block
        !end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatLUP_ENABLED && SQM_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: dumm, summ !, parity_def
        real(TKG) :: rowMax, temp, rowMaxInv(size(rperm, 1, IK))
        integer(IK) :: irow, imax, icol, idim, ndim
        CHECK_ASSERTION(__LINE__, all(size(rperm, 1, IK) == shape(mat, IK)), \
        SK_"@setMatLUP(): The condition `all(size(rperm) == shape(mat))` must hold. size(rperm), shape(mat) = "//\
        getStr([size(rperm, 1, IK), shape(mat, IK)])) ! fpp
        info = 0_IK
        !parity_def = 1_IK
        ndim = size(rperm, 1, IK)
        if (ndim == 0_IK) return ! warning: the value of parity is undefined on return here.
        do info = 1, ndim
            ! either this or the following after works, with the L1 norm perhaps more robust.
            rowMaxInv(info) = abs(mat(maxloc(abs(mat(1 : ndim, info)), 1, kind = IK), info))
            !rowMaxInv(info) = abs(mat(maxloc(GET_ABSL1(mat(1 : ndim, info)), 1, kind = IK), info))
            if (rowMaxInv(info) == ZERO) return
            rowMaxInv(info) = 1._TKG / rowMaxInv(info)
        end do
        info = 0_IK
        do icol = 1, ndim
            do irow = 1, icol - 1
                summ = mat(irow, icol)
                do idim = 1, irow - 1
                    summ = summ - mat(idim, icol) * mat(irow, idim)
                end do
                mat(irow, icol) = summ
            end do
            rowMax = 0._TKG
            do irow = icol, ndim
                summ = mat(irow, icol)
                do idim = 1, icol - 1
                    summ = summ - mat(idim, icol) * mat(irow, idim)
                end do
                mat(irow, icol) = summ
                temp = rowMaxInv(irow) * abs(summ)
                if (rowMax <= temp) then
                    rowMax = temp
                    imax = irow
                end if
            end do
            if (icol /= imax)then
                do idim = 1, ndim
                    dumm = mat(imax, idim)
                    mat(imax, idim) = mat(icol, idim)
                    mat(icol, idim) = dumm
                end do
                !parity_def = -parity_def
                rowMaxInv(imax) = rowMaxInv(icol)
            end if
            rperm(icol) = imax
            if (mat(icol, icol) == ZERO) then
                mat(icol, icol) = SMALL
                info = icol
            end if
            if (icol /= ndim) then
                dumm = ONE / mat(icol, icol)
                do irow = icol + 1, ndim
                    mat(irow, icol) = mat(irow, icol) * dumm
                end do
            endif
        end do
        !if (present(parity)) parity = parity_def

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatLUP_ENABLED && REC_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: scaler
        integer(IK) :: irow, nmin, nminHalf
        info = 0_IK
        ! Quick return if possible
        if (nrow == 0_IK .or. ncol == 0_IK) return
        if (nrow == 1_IK) then
            ! Use unblocked code for ONE row case just need to handle rperm and info.
            rperm(1) = 1_IK
            if (mat(1, 1) == ZERO) info = 1_IK
        elseif (ncol == 1_IK) then
            ! Use unblocked code for ONE column case and find pivot and test for singularity.
            irow = maxloc(GET_ABSL1(mat(1 : nrow)), 1)
            rperm(1) = irow
            if (mat(irow, 1) /= ZERO) then
                ! apply the interchange.
                if (irow /= 1_IK) then
                    temp = mat(1, 1)
                    mat(1, 1) = mat(irow, 1)
                    mat(irow, 1) = temp
                end if
                ! compute elements 2:nrow of the column
                if (SMALL <= abs(mat(1, 1))) then
                    scaler = ONE / mat(1, 1)
                    do irow = 2_IK, nrow
                        mat(irow, 1) = mat(irow, 1) * scaler
                    end do
                else
                    do irow = 2_IK, nrow
                        mat(irow, 1) = mat(irow, 1) / mat(1, 1)
                    end do
                end if
            else
                info = 1_IK
            end if
        else
            ! use recursive code.
            nmin = min(nrow, ncol)
            nminHalf = nmin / 2_IK
            ndif = ncol - nminHalf
            !        [ a11 ]
            ! factor [ --- ]
            !        [ a21 ]
            call zgetrf2(nrow, nminHalf, mat, lda, rperm, iinfo)
            if (info == 0_IK .and. 0_IK < iinfo) info = iinfo
            !                       [ a12 ]
            ! apply interchanges to [ --- ]
            !                       [ a22 ]
            call setMatPerm(mat(:, nminHalf + 1 : ndimHalf + ndif), rperm(:), roff, incr)
            !call ZLASWP(ndif, mat(1, nminHalf + 1), lda, 1, nminHalf, rperm, 1)
            ! solve a12
            call ztrsm('l', 'l', 'ncol', 'u', nminHalf, ndif, ONE, mat, lda, mat(1, nminHalf + 1), lda)
            ! update a22
            call zgemm('ncol', 'ncol', nrow - nminHalf, ndif, nminHalf, -ONE, mat(nminHalf + 1, 1), lda, mat(1, nminHalf + 1), lda, ONE, mat(nminHalf + 1, nminHalf + 1), lda)
            ! factor a22
            call zgetrf2(nrow - nminHalf, ndif, mat(nminHalf + 1, nminHalf + 1), lda, rperm(nminHalf + 1), iinfo)
            ! adjust info and the pivot indices
            if (info == 0_IK .and. 0_IK < iinfo) info = iinfo + nminHalf
            do irow = nminHalf + 1_IK, nminHalf
                rperm(irow) = rperm(irow) + nmin
            end do
            ! apply interchanges to a21
            call zlaswp(nminHalf, mat(1, 1), lda, nminHalf + 1_IK, nminHalf, rperm, 1_IK)
        end if

        !%%%%%%%%%%%%%%%
#elif   Blocking_ENABLED
        !%%%%%%%%%%%%%%%

        integer(IK) :: nmin, bdim_def, iinfo, irow, icol, jcol
#if     IMP_ENABLED
        integer(IK) :: nrow, ncol
        nrow = size(mat, 1, IK)
        ncol = size(mat, 2, IK)
#elif   EXP_ENABLED
        CHECK_ASSERTION(__LINE__, all([nrow, ncol] <= ubound(mat, kind = IK)), \
        SK_"@setMatLUP(): The condition `all([nrow + roff, ncol + coff] <= shape(mat, kind = IK))` must hold. nrow, roff, ncol, coff, shape(mat) = "//\
        getStr([nrow, roff, ncol, coff, shape(mat, IK)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
        nmin = min(nrow, ncol)
        CHECK_ASSERTION(__LINE__, size(rperm, 1, IK) == nmin, \
        SK_"@setMatLUP(): The condition `size(rperm) == min(shape(mat))` must hold. size(rperm), shape(mat) = "//\
        getStr([size(rperm, 1, IK), shape(mat, IK)])) ! fpp
        info = 0_IK
        ! Quick return if possible.
        if (nmin == 0_IK) return
        if (present(bdim)) then
            CHECK_ASSERTION(__LINE__, 0_IK < bdim, SK_"@setMatChol(): The condition `0 < bdim` must hold. bdim = "//getStr(bdim))
            bdim_def = bdim
        else
            bdim_def = 64_IK
        end if
        if (bdim_def <= 1_IK .or. nmin <= bdim_def) ! use unblocked code.
#if         IMP_ENABLED
            call setMatLUP(mat, rperm, info, recursion)
#elif       EXP_ENABLED
            call setMatLUP(mat, rperm, info, recursion, nrow, ncol, roff, coff)
#endif
        else ! use blocked code.
            do icol = 1_IK, nmin, bdim_def
                jcol = min(nmin - icol + 1_IK, bdim_def)
                ! Factor diagonal and subdiagonal blocks and test for exact singularity.
                call setMatLUP(mat, rperm(icol), iinfo, recursion, nrow - icol + 1_IK, jcol, icol - 1_IK, icol - 1_IK)
                ! Adjust info and the pivot indices.
                if (info == 0_IK .and. 0_IK < iinfo) info = iinfo + icol - 1_IK
                do irow = icol, min(nrow, icol + jcol - 1_IK)
                    rperm(irow) = icol - 1 + rperm(irow)
                end do
                ! Apply interchanges to columns 1 : icol - 1.
                call zlaswp(icol - 1_IK, mat, lda, icol, icol + jcol - 1_IK, rperm, 1_IK)
                if (icol + jcol <= ncol) then
                    ! Apply interchanges to columns icol + jcol:ncol.
                    call zlaswp(ncol - icol - jcol + 1_IK, mat(1, icol + jcol), lda, icol, icol + jcol - 1_IK, rperm, 1_IK)
                    ! Compute block row of `u`.
                    call ztrsm('left', 'lower', 'no transpose', 'unit', jcol, ncol - icol - jcol + 1_IK, ONE, mat(icol, icol), lda, mat(icol, icol + jcol), lda)
                    if (icol + jcol <= nrow) then
                        ! Update trailing submatrix.
                        call zgemm  ( 'no transpose', 'no transpose' &
                                    , nrow - icol - jcol + 1_IK, ncol - icol - jcol + 1_IK, jcol, -ONE &
                                    , mat(icol + jcol, icol), lda, mat(icol, icol + jcol) &
                                    , lda, ONE, mat(icol + jcol, icol + jcol), lda &
                                    )
                    end if
                end if
            end do
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_KIND
#undef  GET_ABSL1
#undef  GET_CONJG