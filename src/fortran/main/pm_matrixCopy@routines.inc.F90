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
!>  This include file contains procedure implementation of [pm_matrixCopy](@ref pm_matrixCopy).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the conjugation rule for triangular rfpack formats and Hermitian transpose.
#if     CK_ENABLED && (THO_ENABLED || (TSO_ENABLED && (RFP_RDP_ENABLED || RDP_RFP_ENABLED)))
#define GET_CONJG(X)conjg(X)
#elif   SK_ENABLED || IK_ENABLED || LK_ENABLED || CK_ENABLED || RK_ENABLED
#define GET_CONJG(X)X
#else
#error  "Unrecognized interface."
#endif
!       ! Set the conjugate equation rules for triangular rfpack format.
!#if     RFP_RDP_ENABLED
!#define EQUATE_CONJG(OBJ1, OBJ2, IND1, IND2) OBJ1 IND1 = GET_CONJG(OBJ2 IND2)
!#elif   RDP_RFP_ENABLED
!#define EQUATE_CONJG(OBJ1, OBJ2, IND1, IND2) OBJ1 IND2 = GET_CONJG(OBJ2 IND1)
!#endif
        ! Set the slicing rules for different operations.
#if     AIO_ENABLED
#define GET_SLICE(I, J)I, J
#elif   TSO_ENABLED || THO_ENABLED
#define GET_SLICE(I, J)J, I
#else
#error  "Unrecognized interface."
#endif
        ! Set the runtime check for string length compatibility.
#if     SK_ENABLED
#define CHECK_STRLEN CHECK_ASSERTION(__LINE__, len(source,IK) <= len(destin,IK), \
SK_"@setMatCopy(): The condition `len(source) <= len(destin)` must hold. len(source), len(destin) = "//getStr([len(source,IK), len(destin,IK)])); ! fpp
#else
#define CHECK_STRLEN
#endif
        ! Set the runtime check for `doff` given the subset.
#if     UXD_ENABLED || UXX_ENABLED
#define CHECK_DOFF CHECK_ASSERTION(__LINE__, -nrow <= doff .and. doff <= 0_IK, \
SK_"@setMatCopy(): The condition `-nrow <= doff .and. doff <= 0` must hold. nrow, doff = "//getStr([nrow, doff])); ! fpp
#elif   XLD_ENABLED || XLX_ENABLED
#define CHECK_DOFF CHECK_ASSERTION(__LINE__, 0_IK <= doff .and. doff <= ncol, \
SK_"@setMatCopy(): The condition `0 <= doff .and. doff <= ncol` must hold. doff, ncol = "//getStr([doff, ncol])); ! fpp
#elif   ULX_ENABLED || XXD_ENABLED
#define CHECK_DOFF CHECK_ASSERTION(__LINE__, -nrow <= doff .and. doff <= ncol, \
SK_"@setMatCopy(): The condition `-nrow <= doff .and. doff <= ncol` must hold. nrow, doff, ncol = "//getStr([nrow, doff, ncol])); ! fpp
#else
#error  "Unrecognized interface."
#endif
        ! Define value setting rules for LFP format.
#if     RDP_LFP_ENABLED
#define MATRDP destin
#define MATLFP source
#define EQUATE(OBJ1, OBJ2, IND1, IND2) OBJ1 IND1 = GET_CONJG(OBJ2 IND2)
#elif   LFP_RDP_ENABLED
#define MATLFP destin
#define MATRDP source
#define EQUATE(OBJ1, OBJ2, IND1, IND2) OBJ1 IND2 = GET_CONJG(OBJ2 IND1)
#elif   !(RDP_RDP_ENABLED || RDP_RFP_ENABLED || RFP_RDP_ENABLED)
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getMatCopy_ENABLED && RDP_RDP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(init)) destin = init
#if     AIO_ENABLED
        call setMatCopy(destin, rdpack, source, rdpack, subset, doff)
#elif   TSO_ENABLED || THO_ENABLED
        call setMatCopy(destin, rdpack, source, rdpack, subset, operation, doff)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatCopy_ENABLED && RDP_LFP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: shape_def(2), doff_def
        if (present(doff) .and. present(shape)) then
            shape_def = shape
            doff_def = doff
        else
            CHECK_ASSERTION(__LINE__, .not. (present(shape) .or. present(doff)), \
            SK_"@setMatCopy(): The condition `.not. (present(shape) .or. present(doff))` must hold. present(shape), present(doff) = "\
            //getStr([present(shape), present(doff)])) ! fpp
#if         XXD_ENABLED
            shape_def(1) = size(source, 1, IK)
#elif       UXD_ENABLED || XLD_ENABLED
            shape_def(1) = (getSqrt(size(source, 1, IK) * 8 + 1) - 1) / 2
#endif
            shape_def(2) = shape_def(1)
            doff_def = 0_IK
        end if
        allocate(destin(shape_def(1), shape_def(2)))
        if (present(init)) destin = init
#if     AIO_ENABLED
        call setMatCopy(destin, dpack, source, spack, subset, doff)
#elif   TSO_ENABLED || THO_ENABLED
        call setMatCopy(destin, dpack, source, spack, subset, operation, doff)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatCopy_ENABLED && LFP_RDP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     XXD_ENABLED
        integer(IK) :: doff_def
        doff_def = 0_IK; if (present(doff)) doff_def = doff
        if (doff_def < 0_IK) then
            allocate(destin(min(size(source, 1, IK) + doff_def, size(source, 2, IK))))
        else
            allocate(destin(min(size(source, 1, IK), size(source, 2, IK) - doff_def)))
        end if
#elif   UXD_ENABLED || XLD_ENABLED
        integer(IK) :: dsize
        dsize = 0_IK
        if (all(0_IK < shape(source, IK))) dsize = getMatIndex(dpack, rdpack, shape(source, IK), subset, shape(source, IK), doff)
        allocate(destin(dsize))
#else
#error  "Unrecognized interface."
#endif
        if (present(init)) destin = init
#if     AIO_ENABLED
        call setMatCopy(destin, dpack, source, spack, subset, doff)
#elif   TSO_ENABLED || THO_ENABLED
        call setMatCopy(destin, dpack, source, spack, subset, operation, doff)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatCopy_ENABLED && RDP_RFP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim
        ndim = size(source, 1, IK)
        if (mod(size(source, 2, IK), 2_IK) == 0_IK) ndim = ndim - 1
        allocate(destin(ndim, ndim))
        if (present(init)) destin = init
#if     AIO_ENABLED
        call setMatCopy(destin, dpack, source, spack, subset)
#elif   TSO_ENABLED || THO_ENABLED
#error  "Unrecognized interface."
        !call setMatCopy(destin, dpack, source, spack, subset, operation)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatCopy_ENABLED && RFP_RDP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim
        ndim = size(source, 1, IK)
        if (mod(ndim, 2_IK) == 0_IK) then
            allocate(destin(ndim + 1, ndim / 2))
        else
            allocate(destin(ndim, ndim / 2 + 1))
        end if
        if (present(init)) destin = init
#if     AIO_ENABLED
        call setMatCopy(destin, dpack, source, spack, subset)
#elif   TSO_ENABLED || THO_ENABLED
#error  "Unrecognized interface."
       !call setMatCopy(destin, dpack, source, spack, subset, operation)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatCopy_ENABLED && RDP_RDP_ENABLED && (UXX_ENABLED || XLX_ENABLED || XXD_ENABLED || UXD_ENABLED || XLD_ENABLED || ULX_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: nrow, ncol, icol, doff_def
        nrow = size(source, 1, IK)
        ncol = size(source, 2, IK)
        ! Check the positivity of the lower bounds of `source`.
#if     AIO_ENABLED
        CHECK_ASSERTION(__LINE__, all([nrow, ncol] == shape(destin, IK)), \
        SK_"@setMatCopy(): The condition `all([size(source,1), size(source,2)] == shape(destin, IK))` must hold. shape(source), shape(destin) = "\
        //getStr([nrow, ncol, shape(destin, IK)])) ! fpp
#elif   TSO_ENABLED || THO_ENABLED
        CHECK_ASSERTION(__LINE__, all([ncol, nrow] == shape(destin, IK)), \
        SK_"@setMatCopy(): The condition `all([size(source,2), size(source,1)] == shape(destin, IK))` must hold. shape(source), shape(destin) = "\
        //getStr([nrow, ncol, shape(destin, IK)])) ! fpp
#else
#error  "Unrecognized interface."
#endif
        CHECK_STRLEN
        if (present(doff)) then
            CHECK_DOFF
            doff_def = doff
        else
            doff_def = 0_IK
        end if
#if     UXX_ENABLED || XLX_ENABLED
#define INCREMENT(X,SIGN)X SIGN 1
#elif   XXD_ENABLED || UXD_ENABLED || XLD_ENABLED || ULX_ENABLED
#define INCREMENT(X,SIGN)X
#else
#error  "Unrecognized interface."
#endif
        ! Start the copy action.
        ! The following loops must remain loop due to a macro expansion limitation.
#if     UXD_ENABLED || UXX_ENABLED
        block
            integer(IK) :: limit
            limit = min(nrow + doff_def, ncol)
            do concurrent(icol = 1 : limit)
                destin(GET_SLICE(1 : icol - INCREMENT(doff_def,-), icol)) = GET_CONJG(source(1 : icol - INCREMENT(doff_def,-), icol))
            end do
            do concurrent(icol = limit + 1 : ncol)
                destin(GET_SLICE(1 : nrow, icol)) = GET_CONJG(source(1 : nrow, icol))
            end do
        end block
#elif   XLD_ENABLED || XLX_ENABLED
        do concurrent(icol = 1 : doff_def)
            destin(GET_SLICE(1 : nrow, icol)) = GET_CONJG(source(1 : nrow, icol))
        end do
        do concurrent(icol = 1 : ncol - doff_def)
            destin(GET_SLICE(INCREMENT(icol,+) : nrow, icol + doff_def)) = GET_CONJG(source(INCREMENT(icol,+) : nrow, icol + doff_def))
        end do
#elif   ULX_ENABLED
        block
            integer(IK) :: limit
            if (0_IK == doff_def) then
                limit = min(nrow, ncol)
                do concurrent(icol = 1 : limit)
                    destin(GET_SLICE(1 : icol - 1, icol)) = GET_CONJG(source(1 : icol - 1, icol))
                    destin(GET_SLICE(icol + 1 : nrow, icol)) = GET_CONJG(source(icol + 1 : nrow, icol))
                end do
                do concurrent(icol = limit + 1 : ncol)
                    destin(GET_SLICE(1 : nrow, icol)) = GET_CONJG(source(1 : nrow, icol))
                end do
            elseif (0_IK < doff_def) then
                do concurrent(icol = 1 : doff_def)
                    destin(GET_SLICE(1 : nrow, icol)) = GET_CONJG(source(1 : nrow, icol))
                end do
                do concurrent(icol = 1 : ncol - doff_def)
                    destin(GET_SLICE(1 : icol - 1, icol + doff_def)) = GET_CONJG(source(1 : icol - 1, icol + doff_def))
                    destin(GET_SLICE(icol + 1 : nrow, icol + doff_def)) = GET_CONJG(source(icol + 1 : nrow, icol + doff_def))
                end do
            else
                limit = min(nrow + doff_def, ncol)
                do concurrent(icol = 1 : limit)
                    destin(GET_SLICE(1 : icol - doff_def - 1, icol)) = GET_CONJG(source(1 : icol - doff_def - 1, icol))
                    destin(GET_SLICE(icol - doff_def + 1 : nrow, icol)) = GET_CONJG(source(icol - doff_def + 1 : nrow, icol))
                end do
                do concurrent(icol = limit + 1 : ncol)
                    destin(GET_SLICE(1 : nrow, icol)) = GET_CONJG(source(1 : nrow, icol))
                end do
            end if
        end block
#elif   XXD_ENABLED
        if (0_IK == doff_def) then
            do concurrent(icol = 1 : min(nrow, ncol))
                destin(GET_SLICE(icol, icol)) = GET_CONJG(source(icol, icol))
            end do
        elseif (0_IK < doff_def) then
            do concurrent(icol = 1 : ncol - doff_def)
                destin(GET_SLICE(icol, icol + doff_def)) = GET_CONJG(source(icol, icol + doff_def))
            end do
        else
            do concurrent(icol = 1 : min(nrow + doff_def, ncol))
                destin(GET_SLICE(icol - doff_def, icol)) = GET_CONJG(source(icol - doff_def, icol))
            end do
        end if
#else
#error  "Unrecognized interface."
#endif
#undef  INCREMENT

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatCopy_ENABLED && (RDP_LFP_ENABLED || LFP_RDP_ENABLED) && XXD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: iell
        CHECK_STRLEN
        if (present(doff)) then
            if (doff < 0_IK) then
                CHECK_ASSERTION(__LINE__, size(MATLFP, 1, IK) == min(size(MATRDP, 1, IK) + doff, size(MATRDP, 2, IK)), \
                SK_"@setMatCopy(): The condition `The number of specified diagonal elements in the two matrix subsets must match. doff, shape(destin), shape(source) = "//\
                getStr([doff, shape(destin, IK), shape(source, IK)])) ! fpp
                do concurrent(iell = 1 : size(MATLFP, 1, IK))
                    EQUATE(destin, source, (GET_SLICE(iell - doff, iell)), (iell))
                end do
                return
            elseif (0_IK < doff) then
                CHECK_ASSERTION(__LINE__, \
                size(MATLFP, 1, IK) == min(size(MATRDP, 1, IK), size(MATRDP, 2, IK) - doff), \
                SK_"@setMatCopy(): The condition `The number of specified diagonal elements in the two matrix subsets must match. doff, shape(destin), shape(source) = "//\
                getStr([doff, shape(destin, IK), shape(source, IK)])) ! fpp
                do concurrent(iell = 1 : size(MATLFP, 1, IK))
                    EQUATE(destin, source, (GET_SLICE(iell, iell + doff)), (iell))
                end do
                return
            end if
        else
            CHECK_ASSERTION(__LINE__, size(MATLFP, 1, IK) == minval(shape(MATRDP, IK)), \
            SK_"@setMatCopy(): The condition `The number of specified diagonal elements in the two matrix subsets must match. shape(destin), shape(source) = "//\
            getStr([shape(destin, IK), shape(source, IK)])) ! fpp
            do concurrent(iell = 1 : size(MATLFP, 1, IK))
                EQUATE(destin, source, (GET_SLICE(iell, iell)), (iell))
            end do
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatCopy_ENABLED && (RDP_LFP_ENABLED || LFP_RDP_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: nrow, ncol, irow, icol, nell, doff_def
#if     AIO_ENABLED
        nrow = size(MATRDP, 1, IK)
        ncol = size(MATRDP, 2, IK)
#elif   TSO_ENABLED || THO_ENABLED
        nrow = size(MATRDP, 2, IK)
        ncol = size(MATRDP, 1, IK)
#else
#error  "Unrecognized interface."
#endif
        CHECK_STRLEN
        if (present(doff)) then
            CHECK_DOFF
            doff_def = abs(doff)
        else
            doff_def = 0_IK
        end if
        CHECK_ASSERTION(__LINE__, nrow * ncol - (nrow - doff_def) * (nrow - doff_def - 1) / 2 <= size(MATRDP, kind = IK), \
        SK_"@setMatCopy(): The condition `The number of specified elements in the two matrix subsets must match. doff, shape(destin), shape(source) = "//\
        getStr([doff_def, shape(destin, IK), shape(source, IK)])) ! fpp
        irow = 0_IK
#if     UXD_ENABLED
        do icol = 1, min(nrow - doff_def, ncol)
            nell = icol + doff_def
            EQUATE(destin, source, (GET_SLICE(1 : nell, icol)), (irow + 1 : irow + nell))
            irow = irow + nell
        end do
        do icol = icol, ncol
            nell = irow + nrow
            EQUATE(destin, source, (GET_SLICE(1 : nrow, icol)), (irow + 1 : nell))
            irow = nell
        end do
#elif   XLD_ENABLED && AIO_ENABLED
        do icol = 1, doff_def
            nell = irow + nrow
            EQUATE(destin, source, (1 : nrow, icol), (irow + 1 : nell))
            irow = nell
        end do
        irow = irow + 1
        do icol = 1, ncol - doff_def
            EQUATE(destin, source, (icol : nrow, icol + doff_def), (irow : irow + nrow - icol))
            irow = irow + nrow - icol + 1_IK
        end do
#elif   XLD_ENABLED && (TSO_ENABLED || THO_ENABLED)
        do icol = 1, nrow - doff_def
            nell = doff_def + icol
            EQUATE(destin, source, (icol, 1 : nell), (irow + 1 : irow + nell))
            irow = irow + nell
        end do
        do icol = icol, ncol
            EQUATE(destin, source, (icol, 1 : nrow), (irow + 1 : irow + nrow))
            irow = irow + nrow
        end do
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatCopy_ENABLED && RFP_RDP_ENABLED && AIO_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        logical(LK) :: isNotTransMatB
        integer(IK) :: ndim, ndimHalf, ndimHalfP1 !, irow, icol
        ndim = size(source, 1, IK)

        CHECK_STRLEN

        CHECK_ASSERTION(__LINE__, size(source, 1, IK) == size(source, 2, IK), \
        SK_"@setMatCopy(): The condition `size(source, 1) == size(source, 2)` must hold. shape(source) = "//getStr(shape(source, IK))) ! fpp

        ! quick return if possible.

        if(ndim <= 1_IK) then
            CHECK_ASSERTION(__LINE__, all(shape(destin, IK) == shape(source, IK)), \
            SK_"@setMatCopy(): The condition `all(shape(destin) == shape(source))` must hold. shape(destin), shape(source) = "//getStr([shape(destin, IK), shape(source, IK)])) ! fpp
            if(ndim == 1_IK) destin(1, 1) = GET_CONJG(source(1, 1))
            return
        end if

        ndimHalf = ndim / 2_IK
        ndimHalfP1 = ndimHalf + 1_IK
        isNotTransMatB = logical(size(destin, 2, IK) < size(destin, 1, IK), LK)
        if (ndim < ndimHalf + ndimHalfP1) then ! ndim is even.

            if (isNotTransMatB) then ! `destin` is in normal form (no transposition).
                CHECK_ASSERTION(__LINE__, all(shape(destin, IK) == [size(source, 1, IK) + 1_IK, size(source, 1, IK) / 2_IK]), \
                SK_"@setMatCopy(): The condition `all(shape(destin) == [size(source, 1) + 1, size(source, 1) / 2]))` must hold. shape(destin), shape(source) = "\
                //getStr([shape(destin, IK), shape(source, IK)])) ! fpp
#if             UXD_ENABLED
                call setMatCopy(destin(1 : ndim, 1 : ndimHalf)                , rdpack, source(1 : ndim, ndimHalfP1 : ndim) , rdpack, subset = uppDia, doff = -ndimHalf)
                call setMatCopy(destin(ndimHalf + 2 : ndim + 1, 1 : ndimHalf) , rdpack, source(1 : ndimHalf, 1 : ndimHalf)  , rdpack, subset = uppDia, operation = transHerm)
#elif           XLD_ENABLED
                call setMatCopy(destin(2 : ndim + 1, 1 : ndimHalf), rdpack, source(1 : ndim, 1 : ndimHalf)              , rdpack, subset = lowDia)
                call setMatCopy(destin(1 : ndimHalf, 1 : ndimHalf), rdpack, source(ndimHalfP1 : ndim, ndimHalfP1 : ndim), rdpack, subset = lowDia, operation = transHerm)
#else
#error          "Unrecognized interface."
#endif
            else ! `destin` is in LAPACK-style transposed form. xxx \todo Is the second transposition Hermitian? must be checked with lapack.
                CHECK_ASSERTION(__LINE__, all(shape(destin, IK) == [size(source, 1, IK) / 2_IK, size(source, 1, IK) + 1_IK]), \
                SK_"@setMatCopy(): The condition `all(shape(destin) == [size(source, 1) / 2, size(source, 1) + 1])` must hold. shape(destin), shape(source) = "\
                //getStr([shape(destin, IK), shape(source, IK)])) ! fpp
#if             UXD_ENABLED
                call setMatCopy(destin(1 : ndimHalf, 1 : ndim)                , rdpack, source(1 : ndim, ndimHalfP1 : ndim) , rdpack, subset = uppDia, operation = transHerm, doff = -ndimHalf)
                call setMatCopy(destin(1 : ndimHalf, ndimHalf + 2 : ndim + 1) , rdpack, source(1 : ndimHalf, 1 : ndimHalf)  , rdpack, subset = uppDia)
#elif           XLD_ENABLED
                call setMatCopy(destin(1 : ndimHalf, 2 : ndim + 1), rdpack, source(1 : ndim, 1 : ndimHalf)              , rdpack, subset = lowDia, operation = transHerm)
                call setMatCopy(destin(1 : ndimHalf, 1 : ndimHalf), rdpack, source(ndimHalfP1 : ndim, ndimHalfP1 : ndim), rdpack, subset = lowDia)
#else
#error          "Unrecognized interface."
#endif
            end if

        else ! ndim is odd.

            if (isNotTransMatB) then ! `destin` is in normal form (no transposition).
                CHECK_ASSERTION(__LINE__, all(shape(destin, IK) == [size(source, 1, IK), size(source, 1, IK) / 2 + 1]), \
                SK_"@setMatCopy(): The condition `all(shape(destin) == [size(source, 1), size(source, 1) / 2 + 1]))` must hold. shape(destin), shape(source) = "\
                //getStr([shape(destin, IK), shape(source, IK)])) ! fpp
#if             UXD_ENABLED
                call setMatCopy(destin(1 : ndim, 1 : ndimHalfP1)          , rdpack, source(1 : ndim, ndimHalfP1 : ndim) , rdpack, subset = uppDia, doff = -ndimHalf)
                call setMatCopy(destin(ndimHalf + 2 : ndim, 1 : ndimHalf) , rdpack, source(1 : ndimHalf, 1 : ndimHalf)  , rdpack, subset = uppDia, operation = transHerm)
#elif           XLD_ENABLED
                call setMatCopy(destin(1 : ndim, 1 : ndimHalfP1)      , rdpack, source(1 : ndim, 1 : ndimHalfP1)                    , rdpack, subset = lowDia)
                call setMatCopy(destin(1 : ndimHalf, 2 : ndimHalfP1)  , rdpack, source(ndimHalfP1 + 1 : ndim, ndimHalfP1 + 1 : ndim), rdpack, subset = lowDia, operation = transHerm)
#else
#error          "Unrecognized interface."
#endif
            else ! `destin` is in LAPACK-style transposed form. xxx \todo Is the second transposition Hermitian? must be checked with lapack.
                CHECK_ASSERTION(__LINE__, all(shape(destin, IK) == [size(source, 1, IK) / 2 + 1, size(source, 1, IK)]), \
                SK_"@setMatCopy(): The condition `all(shape(destin) == [size(source, 1) / 2 + 1, size(source, 1)]))` must hold. shape(destin), shape(source) = "\
                //getStr([shape(destin, IK), shape(source, IK)])) ! fpp
#if             UXD_ENABLED
                call setMatCopy(destin(1 : ndimHalfP1, 1 : ndim)          , rdpack, source(1 : ndim, ndimHalfP1 : ndim) , rdpack, subset = uppDia, operation = transHerm, doff = -ndimHalf)
                call setMatCopy(destin(1 : ndimHalf, ndimHalf + 2 : ndim) , rdpack, source(1 : ndimHalf, 1 : ndimHalf)  , rdpack, subset = uppDia)
#elif           XLD_ENABLED
                call setMatCopy(destin(1 : ndimHalfP1, 1 : ndim)      , rdpack, source(1 : ndim, 1 : ndimHalfP1)            , rdpack, subset = lowDia, operation = transHerm)
                call setMatCopy(destin(2 : ndimHalfP1, 1 : ndimHalf)  , rdpack, source(ndimHalfP1 : ndim, ndimHalfP1 : ndim), rdpack, subset = lowDia)
#else
#error          "Unrecognized interface."
#endif
            end if

        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatCopy_ENABLED && RDP_RFP_ENABLED && AIO_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        logical(LK) :: isNotTransMatA
        integer(IK) :: ndim, ndimHalf, ndimHalfP1
        ndim = size(destin, 1, IK)

        CHECK_STRLEN
        CHECK_ASSERTION(__LINE__, size(destin, 1, IK) == size(destin, 2, IK), \
        SK_"@setMatCopy(): The condition `size(destin, 1) == size(destin, 2)` must hold. shape(destin) = "//getStr(shape(destin, IK))) ! fpp

        ! quick return if possible.

        if(ndim <= 1_IK) then
            CHECK_ASSERTION(__LINE__, all(shape(destin, IK) == shape(source, IK)), \
            SK_"@setMatCopy(): The condition `all(shape(destin) == shape(source))` must hold. shape(destin), shape(source) = "//getStr([shape(destin, IK), shape(source, IK)])) ! fpp
            if(ndim == 1_IK) destin(1, 1) = GET_CONJG(source(1, 1))
            return
        end if

        ndimHalf = ndim / 2_IK
        ndimHalfP1 = ndimHalf + 1_IK
        isNotTransMatA = logical(size(source, 2, IK) < size(source, 1, IK), LK)
        if (ndim < ndimHalf + ndimHalfP1) then ! ndim is even.

            if (isNotTransMatA) then ! `source` is in normal form (no transposition).
                CHECK_ASSERTION(__LINE__, all(shape(source, IK) == [size(destin, 1, IK) + 1_IK, size(destin, 1, IK) / 2_IK]), \
                SK_"@setMatCopy(): The condition `all(shape(source) == [size(destin, 1) + 1, size(destin, 1) / 2]))` must hold. shape(source), shape(destin) = "\
                //getStr([shape(source, IK), shape(destin, IK)])) ! fpp
#if             UXD_ENABLED
                call setMatCopy(destin(1 : ndim, ndimHalfP1 : ndim)   , rdpack, source(1 : ndim, 1 : ndimHalf)                  , rdpack, subset = uppDia, doff = -ndimHalf)
                call setMatCopy(destin(1 : ndimHalf, 1 : ndimHalf)    , rdpack, source(ndimHalf + 2 : ndim + 1, 1 : ndimHalf)   , rdpack, subset = uppDia, operation = transHerm)
#elif           XLD_ENABLED
                call setMatCopy(destin(1 : ndim, 1 : ndimHalf)                , rdpack, source(2 : ndim + 1, 1 : ndimHalf), rdpack, subset = lowDia)
                call setMatCopy(destin(ndimHalfP1 : ndim, ndimHalfP1 : ndim)  , rdpack, source(1 : ndimHalf, 1 : ndimHalf), rdpack, subset = lowDia, operation = transHerm)
#else
#error          "Unrecognized interface."
#endif
            else ! `source` is in LAPACK-style transposed form. xxx \todo Is the second transposition Hermitian? must be checked with lapack.
                CHECK_ASSERTION(__LINE__, all(shape(source, IK) == [size(destin, 1, IK) / 2_IK, size(destin, 1, IK) + 1_IK]), \
                SK_"@setMatCopy(): The condition `all(shape(source) == [size(destin, 1) / 2, size(destin, 1) + 1])` must hold. shape(source), shape(destin) = "\
                //getStr([shape(source, IK), shape(destin, IK)])) ! fpp
#if             UXD_ENABLED
                call setMatCopy(destin(1 : ndim, ndimHalfP1 : ndim)   , rdpack, source(1 : ndimHalf, 1 : ndim)                  , rdpack, subset = uppDia, operation = transHerm, doff = -ndimHalf)
                call setMatCopy(destin(1 : ndimHalf, 1 : ndimHalf)    , rdpack, source(1 : ndimHalf, ndimHalf + 2 : ndim + 1)   , rdpack, subset = uppDia)
#elif           XLD_ENABLED
                call setMatCopy(destin(1 : ndim, 1 : ndimHalf)                , rdpack, source(1 : ndimHalf, 2 : ndim + 1), rdpack, subset = lowDia, operation = transHerm)
                call setMatCopy(destin(ndimHalfP1 : ndim, ndimHalfP1 : ndim)  , rdpack, source(1 : ndimHalf, 1 : ndimHalf), rdpack, subset = lowDia)
#else
#error          "Unrecognized interface."
#endif
            end if

        else ! ndim is odd.

            if (isNotTransMatA) then ! `source` is in normal form (no transposition).
                CHECK_ASSERTION(__LINE__, all(shape(source, IK) == [size(destin, 1, IK), size(destin, 1, IK) / 2 + 1]), \
                SK_"@setMatCopy(): The condition `all(shape(source) == [size(destin, 1), size(destin, 1) / 2 + 1]))` must hold. shape(source), shape(destin) = "\
                //getStr([shape(source, IK), shape(destin, IK)])) ! fpp
#if             UXD_ENABLED
                call setMatCopy(destin(1 : ndim, ndimHalfP1 : ndim)   , rdpack, source(1 : ndim, 1 : ndimHalfP1)            , rdpack, subset = uppDia, doff = -ndimHalf)
                call setMatCopy(destin(1 : ndimHalf, 1 : ndimHalf)    , rdpack, source(ndimHalf + 2 : ndim, 1 : ndimHalf)   , rdpack, subset = lowDia, operation = transHerm)
#elif           XLD_ENABLED
                call setMatCopy(destin(1 : ndim, 1 : ndimHalfP1)                      , rdpack, source(1 : ndim, 1 : ndimHalfP1)    , rdpack, subset = lowDia)
                call setMatCopy(destin(ndimHalfP1 + 1 : ndim, ndimHalfP1 + 1 : ndim)  , rdpack, source(1 : ndimHalf, 2 : ndimHalfP1), rdpack, subset = uppDia, operation = transHerm)
#else
#error          "Unrecognized interface."
#endif
            else ! `source` is in LAPACK-style transposed form. xxx \todo Is the second transposition Hermitian? must be checked with lapack.
                CHECK_ASSERTION(__LINE__, all(shape(source, IK) == [size(destin, 1, IK) / 2 + 1, size(destin, 1, IK)]), \
                SK_"@setMatCopy(): The condition `all(shape(source) == [size(destin, 1) / 2 + 1, size(destin, 1)]))` must hold. shape(source), shape(destin) = "\
                //getStr([shape(source, IK), shape(destin, IK)])) ! fpp
#if             UXD_ENABLED
                call setMatCopy(destin(1 : ndim, ndimHalfP1 : ndim)   , rdpack, source(1 : ndimHalfP1, 1 : ndim)            , rdpack, subset = uppDia, operation = transHerm, doff = -ndimHalf)
                call setMatCopy(destin(1 : ndimHalf, 1 : ndimHalf)    , rdpack, source(1 : ndimHalf, ndimHalf + 2 : ndim)   , rdpack, subset = uppDia)
#elif           XLD_ENABLED
                call setMatCopy(destin(1 : ndim, 1 : ndimHalfP1)              , rdpack, source(1 : ndimHalfP1, 1 : ndim)    , rdpack, subset = lowDia, operation = transHerm)
                call setMatCopy(destin(ndimHalfP1 : ndim, ndimHalfP1 : ndim)  , rdpack, source(2 : ndimHalfP1, 1 : ndimHalf), rdpack, subset = lowDia)
#else
#error          "Unrecognized interface."
#endif
            end if

        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatCopy_ENABLED && RFP_RDP_ENABLED && (TSO_ENABLED || THO_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#error  "Unrecognized interface."

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatCopy_ENABLED && RDP_RFP_ENABLED && (TSO_ENABLED || THO_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#error  "Unrecognized interface."

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef CHECK_STRLEN
#undef CHECK_DOFF
#undef GET_CONJG
#undef GET_SLICE
#undef MATLFP
#undef MATRDP
#undef EQUATE