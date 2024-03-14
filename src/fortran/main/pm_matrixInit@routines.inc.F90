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
!>  This include file contains procedure implementation of the generic interface [pm_matrixInit](@ref pm_matrixInit).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Thursday 01:00 AM, September 23, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getMatInit_ENABLED && XXD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: roff_def, coff_def

        ! set the row offset.
        if (present(roff)) then
            roff_def = roff
        else
            roff_def = 0_IK
        end if

        ! set the column offset.
        if (present(coff)) then
            coff_def = coff
        else
            coff_def = 0_IK
        end if

#if     D2XX0_ENABLED
        ! set the number of the diagonal elements initialized.
        call setMatInit(mat, subset, vdia, getOption(min(shape(1) - roff_def, shape(2) - coff_def), ndia), roff_def, coff_def)
#elif   D2XX1_ENABLED
        call setMatInit(mat, subset, vdia, roff_def, coff_def)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatInit_ENABLED && XLD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: roff_def, coff_def
        integer(IK) :: nrow_def, ncol_def

        ! set the row offset.
        if (present(roff)) then
            roff_def = roff
        else
            roff_def = 0_IK
        end if

        ! set the column offset.
        if (present(coff)) then
            coff_def = coff
        else
            coff_def = 0_IK
        end if

        ! set the number of the rows initialized.
        if (present(nrow)) then
            nrow_def = nrow
        else
            nrow_def = shape(1) - roff_def
        end if

        ! set the number of the columns initialized.
        if (present(ncol)) then
            ncol_def = ncol
        else
            ncol_def = shape(2) - coff_def
        end if

        ! set the number of the columns initialized.
        if (present(ncol)) then
            ncol_def = ncol
        else
            ncol_def = shape(2) - coff_def
        end if

        call setMatInit(mat, subset, vlow, vdia, nrow_def, ncol_def, roff_def, coff_def, doff)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatInit_ENABLED && UXD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: roff_def, coff_def
        integer(IK) :: nrow_def, ncol_def

        ! set the row offset.
        if (present(roff)) then
            roff_def = roff
        else
            roff_def = 0_IK
        end if

        ! set the column offset.
        if (present(coff)) then
            coff_def = coff
        else
            coff_def = 0_IK
        end if

        ! set the number of the rows initialized.
        if (present(nrow)) then
            nrow_def = nrow
        else
            nrow_def = shape(1) - roff_def
        end if

        ! set the number of the columns initialized.
        if (present(ncol)) then
            ncol_def = ncol
        else
            ncol_def = shape(2) - coff_def
        end if

        ! set the number of the columns initialized.
        if (present(ncol)) then
            ncol_def = ncol
        else
            ncol_def = shape(2) - coff_def
        end if

        call setMatInit(mat, subset, vupp, vdia, nrow_def, ncol_def, roff_def, coff_def, doff)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatInit_ENABLED && ULX_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: roff_def, coff_def
        integer(IK) :: nrow_def, ncol_def

        ! set the row offset.
        if (present(roff)) then
            roff_def = roff
        else
            roff_def = 0_IK
        end if

        ! set the column offset.
        if (present(coff)) then
            coff_def = coff
        else
            coff_def = 0_IK
        end if

        ! set the number of the rows initialized.
        if (present(nrow)) then
            nrow_def = nrow
        else
            nrow_def = shape(1) - roff_def
        end if

        ! set the number of the columns initialized.
        if (present(ncol)) then
            ncol_def = ncol
        else
            ncol_def = shape(2) - coff_def
        end if

        ! set the number of the columns initialized.
        if (present(ncol)) then
            ncol_def = ncol
        else
            ncol_def = shape(2) - coff_def
        end if

        call setMatInit(mat, subset, vupp, vlow, nrow_def, ncol_def, roff_def, coff_def, doff)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMatInit_ENABLED && ULD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: roff_def, coff_def
        integer(IK) :: nrow_def, ncol_def

        ! set the row offset.
        if (present(roff)) then
            roff_def = roff
        else
            roff_def = 0_IK
        end if

        ! set the column offset.
        if (present(coff)) then
            coff_def = coff
        else
            coff_def = 0_IK
        end if

        ! set the number of the rows initialized.
        if (present(nrow)) then
            nrow_def = nrow
        else
            nrow_def = shape(1) - roff_def
        end if

        ! set the number of the columns initialized.
        if (present(ncol)) then
            ncol_def = ncol
        else
            ncol_def = shape(2) - coff_def
        end if

        ! set the number of the columns initialized.
        if (present(ncol)) then
            ncol_def = ncol
        else
            ncol_def = shape(2) - coff_def
        end if

        call setMatInit(mat, subset, vupp, vlow, vdia, nrow_def, ncol_def, roff_def, coff_def, doff)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatInit_ENABLED && XXD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
#if     IMP_ENABLED
        integer(IK), parameter :: roff = 0_IK, coff = 0_IK
#elif   !EXP_ENABLED
#error  "Unrecognized interface."
#endif
#if     D2XXF_ENABLED
        integer(IK) :: ndia
        ndia = min(size(mat, 1, IK) - roff, size(mat, 2, IK) - coff)
#define GET_DIA(i) vdia
#elif   D2XX0_ENABLED
#define NDIA ndia
#define GET_DIA(i) vdia
        CHECK_ASSERTION(__LINE__, ndia <= minval(shape(mat, IK) - [roff, coff]), SK_"@setMatInit(): The condition `ndia <= minval(shape(mat, IK) - [roff, coff])` must hold. ndia, shape(mat), roff, coff = "//getStr([ndia, shape(mat, kind = IK), roff, coff]))
#elif   D2XX1_ENABLED
#define NDIA size(vdia, kind = IK)
#define GET_DIA(i) vdia(i)
        CHECK_ASSERTION(__LINE__, NDIA <= minval(shape(mat, IK) - [roff, coff]), SK_"@setMatInit(): The condition `size(vdia) <= minval(shape(mat, IK) - [roff, coff])` must hold. size(vdia), shape(mat), roff, coff = "//getStr([NDIA, shape(mat, kind = IK), roff, coff]))
#else
#error  "Unrecognized interface."
#endif
#if     SK_ENABLED
        CHECK_ASSERTION(__LINE__, len(vdia, IK) <= len(mat, IK), \
        SK_"@setMatInit(): The condition `len(vdia) <= len(mat)` must hold. len(vdia), len(mat) = "//getStr([len(vdia, IK), len(mat, IK)]))
#endif
        CHECK_ASSERTION(__LINE__, 0_IK <= roff, SK_"@setMatInit(): The condition `0_IK <= roff` must hold. roff = "//getStr(roff))
        CHECK_ASSERTION(__LINE__, 0_IK <= coff, SK_"@setMatInit(): The condition `0_IK <= coff` must hold. coff = "//getStr(coff))
        do concurrent(i = 1 : NDIA)
            mat(i + roff, i + coff) = GET_DIA(i)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatInit_ENABLED && (XLX_ENABLED || XLD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenDia, doff_def
#if     IMP_ENABLED
        integer(IK), parameter :: roff = 0_IK, coff = 0_IK
        integer(IK) :: nrow, ncol
        nrow = size(mat, 1, IK)
        ncol = size(mat, 2, IK)
#elif   EXP_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK <= roff, SK_"@setMatInit(): The condition `0_IK <= roff` must hold. roff = "//getStr(roff))
        CHECK_ASSERTION(__LINE__, 0_IK <= coff, SK_"@setMatInit(): The condition `0_IK <= coff` must hold. coff = "//getStr(coff))
        CHECK_ASSERTION(__LINE__, 0_IK <= nrow .and. nrow + roff <= size(mat, 1, IK), SK_"@setMatInit(): The condition `0_IK <= nrow .and. nrow + roff <= size(mat,1)` must hold. nrow, roff, size(mat,1) = "//getStr([nrow, roff, size(mat, 1, IK)]))
        CHECK_ASSERTION(__LINE__, 0_IK <= ncol .and. ncol + coff <= size(mat, 2, IK), SK_"@setMatInit(): The condition `0_IK <= ncol .and. ncol + coff <= size(mat,2)` must hold. ncol, coff, size(mat,2) = "//getStr([ncol, coff, size(mat, 2, IK)]))
#else
#error  "Unrecognized interface."
#endif
#if     SK_ENABLED
        CHECK_ASSERTION(__LINE__, len(vlow, IK) <= len(mat, IK), SK_"@setMatInit(): The condition `len(vlow) <= len(mat)` must hold. len(vlow), len(mat) = "//getStr([len(vlow, IK), len(mat, IK)]))
#if     XLD_ENABLED
        CHECK_ASSERTION(__LINE__, len(vdia, IK) <= len(mat, IK), SK_"@setMatInit(): The condition `len(vdia) <= len(mat)` must hold. len(vdia), len(mat) = "//getStr([len(vdia, IK), len(mat, IK)]))
#endif
#endif
        if (nrow < 1_IK .or. ncol < 1_IK) return ! this must be here.
        if (present(doff)) then
            CHECK_ASSERTION(__LINE__, 0 <= doff .and. doff < ncol, SK_"@setMatInit(): The condition `0 <= doff .and. doff < ncol` must hold. doff, ncol = "//getStr([doff, ncol]))
            doff_def = doff
        else
            doff_def = 0_IK
        endif
        CHECK_ASSERTION(__LINE__, ncol - doff_def <= nrow, SK_"@setMatInit(): The condition `ncol - doff <= nrow` must hold. You cannot ask for more columns than possible. ncol, doff, nrow = "//getStr([ncol, doff_def, nrow]))
#if     D2X00_ENABLED || D2X0X_ENABLED
#define GET_DIA(i) vdia
        lenDia = ncol - doff_def ! min(ncol - doff_def, nrow)
#elif   D2X01_ENABLED
#define GET_DIA(i) vdia(i)
        lenDia = size(vdia, 1, IK)
        CHECK_ASSERTION(__LINE__, lenDia == ncol - doff_def, SK_"@setMatInit(): The condition `size(vdia) == ncol - doff` must hold. size(vdia), nrow, doff, ncol = "//getStr([lenDia, nrow, doff_def, ncol]))
#else
#error  "Unrecognized interface."
#endif
        !check_assertion(__LINE__, all(0_IK < [nrow, ncol]), SK_"@setMatInit(): The condition `all(0 < [nrow, ncol])` must hold. nrow, ncol = "//getStr([[nrow, ncol]]))
        do concurrent(i = 1 : doff_def)
            mat(1 : nrow, i) = vlow
        end do
#if     XLD_ENABLED
        mat(1 , 1 + doff_def) = GET_DIA(1)
#endif
        mat(2 : nrow , 1 + doff_def) = vlow
        do i = 2_IK, lenDia
#if         XLD_ENABLED
            mat(i , i + doff_def) = GET_DIA(i)
#endif
            mat(i + 1 : nrow , i + doff_def) = vlow
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatInit_ENABLED && (UXX_ENABLED || UXD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenDia, doff_def
#if     IMP_ENABLED
        integer(IK), parameter :: roff = 0_IK, coff = 0_IK
        integer(IK) :: nrow, ncol
        nrow = size(mat, 1, IK)
        ncol = size(mat, 2, IK)
#elif   EXP_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK <= roff, SK_"@setMatInit(): The condition `0_IK <= roff` must hold. roff = "//getStr(roff))
        CHECK_ASSERTION(__LINE__, 0_IK <= coff, SK_"@setMatInit(): The condition `0_IK <= coff` must hold. coff = "//getStr(coff))
        CHECK_ASSERTION(__LINE__, 0_IK <= nrow .and. nrow + roff <= size(mat, 1, IK), SK_"@setMatInit(): The condition `0_IK <= nrow .and. nrow + roff <= size(mat,1)` must hold. nrow, roff, size(mat,1) = "//getStr([nrow, roff, size(mat, 1, IK)]))
        CHECK_ASSERTION(__LINE__, 0_IK <= ncol .and. ncol + coff <= size(mat, 2, IK), SK_"@setMatInit(): The condition `0_IK <= ncol .and. ncol + coff <= size(mat,2)` must hold. ncol, coff, size(mat,2) = "//getStr([ncol, coff, size(mat, 2, IK)]))
#else
#error  "Unrecognized interface."
#endif
#if     SK_ENABLED
        CHECK_ASSERTION(__LINE__, len(vupp, IK) <= len(mat, IK), SK_"@setMatInit(): The condition `len(vupp) <= len(mat)` must hold. len(vupp), len(mat) = "//getStr([len(vupp, IK), len(mat, IK)]))
#if     UXD_ENABLED
        CHECK_ASSERTION(__LINE__, len(vdia, IK) <= len(mat, IK), SK_"@setMatInit(): The condition `len(vdia) <= len(mat)` must hold. len(vdia), len(mat) = "//getStr([len(vdia, IK), len(mat, IK)]))
#endif
#endif
        !check_assertion(__LINE__, all(0_IK < [nrow, ncol]), SK_"@setMatInit(): The condition `all(0 < [nrow, ncol])` must hold. nrow, ncol = "//getStr([[nrow, ncol]]))
        if (nrow < 1_IK .or. ncol < 1_IK) return ! this must be here.
        if (present(doff)) then
            CHECK_ASSERTION(__LINE__, -nrow < doff .and. doff <= 0_IK, SK_"@setMatInit(): The condition `-nrow < doff .and. doff <= 0_IK` must hold. nrow, doff = "//getStr([nrow, doff]))
            doff_def = -doff
        else
            doff_def = 0_IK
        endif
        CHECK_ASSERTION(__LINE__, nrow - doff_def <= ncol, SK_"@setMatInit(): The condition `nrow + doff <= ncol` must hold. You cannot ask for more rows than possible. nrow, doff, ncol = "//getStr([nrow, -doff_def, ncol]))
#if     D20X0_ENABLED || D20XX_ENABLED
#define GET_DIA(i) vdia
        lenDia = nrow - doff_def ! min(nrow - doff_def, ncol)
#elif   D20X1_ENABLED
        lenDia = size(vdia, 1, IK)
        CHECK_ASSERTION(__LINE__, lenDia == nrow - doff_def, SK_"@setMatInit(): The condition `size(vdia) == nrow + doff` must hold. size(vdia), nrow, doff = "//getStr([lenDia, nrow, doff_def]))
#define GET_DIA(i) vdia(i)
#else
#error  "Unrecognized interface."
#endif
        do i = 1_IK, lenDia
            mat(1 : i + doff_def - 1, i) = vupp
#if         UXD_ENABLED
            mat(i + doff_def , i) = GET_DIA(i)
#endif
        end do
        do i = lenDia + 1, ncol
            mat(1 : nrow, i) = vupp
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatInit_ENABLED && ULX_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenDia, doff_def
#if     IMP_ENABLED
        integer(IK), parameter :: roff = 0_IK, coff = 0_IK
        integer(IK) :: nrow, ncol
        nrow = size(mat, 1, IK)
        ncol = size(mat, 2, IK)
#elif   EXP_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK <= roff, SK_"@setMatInit(): The condition `0_IK <= roff` must hold. roff = "//getStr(roff))
        CHECK_ASSERTION(__LINE__, 0_IK <= coff, SK_"@setMatInit(): The condition `0_IK <= coff` must hold. coff = "//getStr(coff))
        CHECK_ASSERTION(__LINE__, 0_IK <= nrow .and. nrow + roff <= size(mat, 1, IK), SK_"@setMatInit(): The condition `0_IK <= nrow .and. nrow + roff <= size(mat,1)` must hold. nrow, roff, size(mat,1) = "//getStr([nrow, roff, size(mat, 1, IK)]))
        CHECK_ASSERTION(__LINE__, 0_IK <= ncol .and. ncol + coff <= size(mat, 2, IK), SK_"@setMatInit(): The condition `0_IK <= ncol .and. ncol + coff <= size(mat,2)` must hold. ncol, coff, size(mat,2) = "//getStr([ncol, coff, size(mat, 2, IK)]))
#else
#error  "Unrecognized interface."
#endif
#if     SK_ENABLED
        CHECK_ASSERTION(__LINE__, len(vupp, IK) <= len(mat, IK), SK_"@setMatInit(): The condition `len(vupp) <= len(mat)` must hold. len(vupp), len(mat) = "//getStr([len(vupp, IK), len(mat, IK)]))
        CHECK_ASSERTION(__LINE__, len(vlow, IK) <= len(mat, IK), SK_"@setMatInit(): The condition `len(vlow) <= len(mat)` must hold. len(vlow), len(mat) = "//getStr([len(vlow, IK), len(mat, IK)]))
#endif
        !check_assertion(__LINE__, all(0_IK < [nrow, ncol]), SK_"@setMatInit(): The condition `all(0 < [nrow, ncol])` must hold. nrow, ncol = "//getStr([[nrow, ncol]]))
        if (nrow < 1_IK .or. ncol < 1_IK) return ! this must be here.
        if (present(doff)) then
            CHECK_ASSERTION(__LINE__, -nrow < doff .and. doff < ncol, SK_"@setMatInit(): The condition `-nrow < doff .and. doff < ncol` must hold. nrow, doff, ncol = "//getStr([nrow, doff, ncol]))
            doff_def = doff
            if (doff_def < 0_IK) then
                lenDia = min(nrow + doff_def, ncol)
            else
                lenDia = min(ncol - doff_def, nrow)
            end if
        else
            doff_def = 0_IK
            lenDia = min(ncol, nrow)
        endif
        if (doff_def > 0_IK) then
            do i = 1, doff_def
                mat(1 : nrow, i) = vlow
            end do
            mat(2 : nrow , doff_def + 1_IK) = vlow
            do i = 2_IK, lenDia
                mat(1 : i - 1, i + doff_def) = vupp
                mat(i + 1 : nrow , i + doff_def) = vlow
            end do
            do i = lenDia + doff_def + 1, ncol
                mat(1 : nrow, i) = vupp
            end do
        else
            do i = 1_IK, lenDia
                mat(1 : i - doff_def - 1, i) = vupp
                mat(i - doff_def + 1 : nrow , i) = vlow
            end do
            do i = lenDia + 1, ncol
                mat(1 : nrow, i) = vupp
            end do
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setMatInit_ENABLED && ULD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenDia, doff_def
#if     IMP_ENABLED
        integer(IK), parameter :: roff = 0_IK, coff = 0_IK
        integer(IK) :: nrow, ncol
        nrow = size(mat, 1, IK)
        ncol = size(mat, 2, IK)
#elif   EXP_ENABLED
        CHECK_ASSERTION(__LINE__, 0_IK <= roff, SK_"@setMatInit(): The condition `0_IK <= roff` must hold. roff = "//getStr(roff))
        CHECK_ASSERTION(__LINE__, 0_IK <= coff, SK_"@setMatInit(): The condition `0_IK <= coff` must hold. coff = "//getStr(coff))
        CHECK_ASSERTION(__LINE__, 0_IK <= nrow .and. nrow + roff <= size(mat, 1, IK), SK_"@setMatInit(): The condition `0_IK <= nrow .and. nrow + roff <= size(mat,1)` must hold. nrow, roff, size(mat,1) = "//getStr([nrow, roff, size(mat, 1, IK)]))
        CHECK_ASSERTION(__LINE__, 0_IK <= ncol .and. ncol + coff <= size(mat, 2, IK), SK_"@setMatInit(): The condition `0_IK <= ncol .and. ncol + coff <= size(mat,2)` must hold. ncol, coff, size(mat,2) = "//getStr([ncol, coff, size(mat, 2, IK)]))
#else
#error  "Unrecognized interface."
#endif
#if     SK_ENABLED
        CHECK_ASSERTION(__LINE__, len(vupp, IK) <= len(mat, IK), SK_"@setMatInit(): The condition `len(vupp) <= len(mat)` must hold. len(vupp), len(mat) = "//getStr([len(vupp, IK), len(mat, IK)]))
        CHECK_ASSERTION(__LINE__, len(vlow, IK) <= len(mat, IK), SK_"@setMatInit(): The condition `len(vlow) <= len(mat)` must hold. len(vlow), len(mat) = "//getStr([len(vlow, IK), len(mat, IK)]))
        CHECK_ASSERTION(__LINE__, len(vdia, IK) <= len(mat, IK), SK_"@setMatInit(): The condition `len(vdia) <= len(mat)` must hold. len(vdia), len(mat) = "//getStr([len(vdia, IK), len(mat, IK)]))
#endif
        !check_assertion(__LINE__, all(0_IK < [nrow, ncol]), SK_"@setMatInit(): The condition `all(0 < [nrow, ncol])` must hold. nrow, ncol = "//getStr([[nrow, ncol]]))
        if (nrow < 1_IK .or. ncol < 1_IK) return ! this must be here.
        if (present(doff)) then
            CHECK_ASSERTION(__LINE__, -nrow < doff .and. doff < ncol, SK_"@setMatInit(): The condition `-nrow < doff .and. doff < ncol` must hold. nrow, doff, ncol = "//getStr([nrow, doff, ncol]))
            doff_def = doff
        else
            doff_def = 0_IK
        endif
#if     D2000_ENABLED
#define GET_DIA(i) vdia
        if (doff_def < 0_IK) then
            lenDia = min(nrow + doff_def, ncol)
        else
            lenDia = min(ncol - doff_def, nrow)
        end if
#elif   D2001_ENABLED
#define GET_DIA(i) vdia(i)
        lenDia = size(vdia, 1, IK)
        CHECK_ASSERTION(__LINE__, lenDia == merge(min(nrow + doff_def, ncol), min(nrow, ncol - doff_def), doff_def < 0_IK), \
        SK_"@setMatInit(): The condition `size(vdia) == merge(min(nrow + doff, ncol), min(nrow, ncol - doff), doff < 0)` must hold. size(vdia), nrow, ncol, doff = "//\
        getStr([size(vdia, 1, IK), nrow, ncol, doff_def]))
#else
#error  "Unrecognized interface."
#endif
        !check_assertion(__LINE__, all(0_IK < [nrow, ncol]), SK_"@setMatInit(): The condition `all(0 < [nrow, ncol])` must hold. nrow, ncol = "//getStr([[nrow, ncol]]))
        if (nrow < 1_IK .or. ncol < 1_IK) return ! this must be here.
        if (doff_def > 0_IK) then
            do i = 1, doff_def
                mat(1 : nrow, i) = vlow
            end do
            i = doff_def + 1_IK
            mat(1 , i) = GET_DIA(1)
            mat(2 : nrow , i) = vlow
            do i = 2_IK, lenDia
                mat(1 : i - 1, i + doff_def) = vupp
                mat(i , i + doff_def) = GET_DIA(i)
                mat(i + 1 : nrow , i + doff_def) = vlow
            end do
            do i = lenDia + doff_def + 1, ncol
                mat(1 : nrow, i) = vupp
            end do
        else
            do i = 1_IK, lenDia
                mat(1 : i - doff_def - 1, i) = vupp
                mat(i - doff_def , i) = GET_DIA(i)
                mat(i - doff_def + 1 : nrow , i) = vlow
            end do
            do i = lenDia + 1, ncol
                mat(1 : nrow, i) = vupp
            end do
        end if
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_DIA
#undef  NDIA