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
!>  This file contains procedure implementations of [pm_matrixIndex](@ref pm_matrixIndex).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !getMatIndex_ENABLED
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     AIO_ENABLED && (UXD_ENABLED || XLD_ENABLED) && (LFP_RDP_ENABLED || RDP_LFP_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim, doffAbs
        CHECK_ASSERTION(__LINE__, all(0_IK <= shape), \
        SK_"@getMatIndex(): The condition `all(0_IK <= shape)` must hold. shape, doff = "//getStr(shape))
        CHECK_ASSERTION(__LINE__, all([0_IK < sindex]), SK_"@getMatIndex(): The condition `all([0 < sindex])` must hold. sindex = "//getStr(sindex))
        if (present(doff)) then
#if         UXD_ENABLED
            doffAbs = -doff
            CHECK_ASSERTION(__LINE__, 0_IK <= shape(1) + doff .and. doff <= 0_IK, \
            SK_"@getMatIndex(): The condition `0 <= shape(1) + doff .and. doff <= 0` must hold. shape, doff = "//getStr([shape, doff]))
#elif       XLD_ENABLED
            doffAbs = doff
            CHECK_ASSERTION(__LINE__, 0_IK <= shape(2) - doff .and. 0_IK <= doff, \
            SK_"@getMatIndex(): The condition `0 <= shape(2) - doff .and. 0 <= doff` must hold. shape, doff = "//getStr([shape, doff]))
#else
#error      "Unrecognized interface."
#endif
        else
            doffAbs = 0_IK
        end if
#if     LFP_RDP_ENABLED && UXD_ENABLED
        ndim = min(sindex(2), shape(1) - doffAbs) ! the effective triangle rank.
        CHECK_ASSERTION(__LINE__, sindex(1) <= sindex(2) + doffAbs, \
        SK_"@getMatIndex(): The condition `sindex(1) <= sindex(2) - doff` must hold. sindex, doff = "//getStr([sindex, -doffAbs]))
        CHECK_ASSERTION(__LINE__, all(sindex <= shape), \
        SK_"@getMatIndex(): The condition `all(sindex <= shape)` must hold. sindex, shape = "//getStr([sindex, shape]))
        dindex & ! LCOV_EXCL_LINE
        = sindex(2) * doffAbs & ! The top rectangle ! LCOV_EXCL_LINE
        + ndim * (ndim - 1_IK) / 2_IK & ! the bottom upper triangle ! LCOV_EXCL_LINE
        + (sindex(2) - ndim) * ndim & ! The rightmost rectangle. ! LCOV_EXCL_LINE
        + sindex(1) - doffAbs ! last column.
#elif   LFP_RDP_ENABLED && XLD_ENABLED
        ndim = shape(2) - doffAbs ! empty triangle rank.
        CHECK_ASSERTION(__LINE__, sindex(2) <= sindex(1) + doffAbs, \
        SK_"@getMatIndex(): The condition `sindex(2) <= sindex(1) + doff` must hold. sindex, doff = "//getStr([sindex, doffAbs]))
        CHECK_ASSERTION(__LINE__, all(sindex <= shape), \
        SK_"@getMatIndex(): The condition `all(sindex <= shape)` must hold. sindex, shape = "//getStr([sindex, shape]))
        block
            integer(IK) :: jcol
            jcol = max(0_IK, sindex(2) - doffAbs - 1_IK)
            dindex & ! LCOV_EXCL_LINE
            = shape(1) * max(0_IK, sindex(2) - 1_IK) & ! The leftmost full rectangle. ! LCOV_EXCL_LINE
            - jcol * (jcol - 1_IK) / 2_IK & ! The empty upper triangle. ! LCOV_EXCL_LINE
            + sindex(1) - jcol
        end block
#elif   RDP_LFP_ENABLED && UXD_ENABLED
#define NDIM min(shape(2), shape(1) - doffAbs)
        CHECK_ASSERTION(__LINE__, sindex <= product(shape) - NDIM * (NDIM - 1_IK) / 2_IK, \
        SK_"@getMatIndex(): The condition `sindex <= product(shape) - NDIM * (NDIM - 1) / 2` must hold. sindex, shape, NDIM = "\
        //getStr([sindex, shape, NDIM]))
#undef  NDIM
        dindex(1) = sindex
        dindex(2) = 1_IK
        do
            ndim = min(shape(1), dindex(2) + doffAbs)
            if (dindex(1) - ndim <= 0_IK) exit
            dindex(1) = dindex(1) - ndim
            dindex(2) = dindex(2) + 1_IK
        end do
#elif   RDP_LFP_ENABLED && XLD_ENABLED
        ndim = shape(2) - doffAbs
        CHECK_ASSERTION(__LINE__, sindex <= product(shape) - ndim * (ndim - 1_IK) / 2_IK, \
        SK_"@getMatIndex(): The condition `sindex <= product(shape) - ndim * (ndim - 1) / 2` must hold. sindex, shape, ndim = "\
        //getStr([sindex, shape, ndim]))
        dindex(1) = sindex
        dindex(2) = 1_IK
        do
            if (doffAbs < dindex(2)) exit
            if (dindex(1) - shape(1) <= 0_IK) exit
            dindex(1) = dindex(1) - shape(1)
            dindex(2) = dindex(2) + 1_IK
        end do
        ndim = shape(1)
        do
            if (dindex(1) - ndim <= 0_IK) then
                dindex(1) = dindex(1) + shape(1) - ndim
                exit
            end if
            dindex(1) = dindex(1) - ndim
            dindex(2) = dindex(2) + 1_IK
            ndim = ndim - 1_IK
            CHECK_ASSERTION(__LINE__, \
            0_IK <= ndim, SK_"@getMatIndex(): Internal error occurred. The condition `0 <= ndim` must hold. ndim = "//getStr(ndim))
        end do
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   AIO_ENABLED && UXD_ENABLED && RFP_RDP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim, ndimHalf
        CHECK_ASSERTION(__LINE__, all(0_IK <= shape), \
        SK_"@getMatIndex(): The condition `all(0_IK <= shape)` must hold. shape, doff = "//getStr(shape))
        ndim = shape(1)
        ndimHalf = ndim / 2_IK
        CHECK_ASSERTION(__LINE__, shape(1) == shape(2), \
        SK_"@getMatIndex(): The condition `shape(1) == shape(2)` must hold. shape = "//getStr(shape))
        CHECK_ASSERTION(__LINE__, 1_IK <= sindex(1) .and. sindex(1) <= sindex(2), \
        SK_"@getMatIndex(): The condition `1 <= sindex(1) .and. sindex(1) <= sindex(2)` must hold. sindex = "//getStr(sindex))
        CHECK_ASSERTION(__LINE__, all(0_IK < sindex) .and. all(sindex <= shape), \
        SK_"@getMatIndex(): The condition `all(0_IK < sindex) .and. all(sindex <= shape)` must hold. sindex, shape = "//getStr([sindex, shape]))
        if (sindex(2) <= ndimHalf) then
            dindex(1) = sindex(2) + ndimHalf + 1_IK
            dindex(2) = sindex(1)
        else
            dindex(1) = sindex(1)
            dindex(2) = sindex(2) - ndimHalf
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   AIO_ENABLED && UXD_ENABLED && RDP_RFP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim, ndimHalf, remainder
        CHECK_ASSERTION(__LINE__, all(0_IK <= shape), \
        SK_"@getMatIndex(): The condition `all(0_IK <= shape)` must hold. shape, doff = "//getStr(shape))
        ndim = shape(1)
        ndimHalf = ndim / 2_IK
        remainder = ndim - ndimHalf * 2_IK
        CHECK_ASSERTION(__LINE__, shape(1) == shape(2), \
        SK_"@getMatIndex(): The condition `shape(1) == shape(2)` must hold. shape = "//getStr(shape))
        CHECK_ASSERTION(__LINE__, 1_IK <= sindex(2) .and. sindex(2) <= ndim + 1_IK - remainder, \
        SK_"@getMatIndex(): The condition `1 <= sindex(1) .and. sindex(1) <= ndim + 1 - mod(ndim, 2)` must hold. sindex = "//getStr(sindex))
        CHECK_ASSERTION(__LINE__, 1_IK <= sindex(1) .and. sindex(2) <= ndimHalf + remainder, \
        SK_"@getMatIndex(): The condition `1 <= sindex(2) .and. sindex(2) <= ndimHalf + mod(ndim, 2)` must hold. sindex = "//getStr(sindex))
        if (sindex(1) - sindex(2) <= ndimHalf) then
            dindex(1) = sindex(1)
            dindex(2) = sindex(2) + ndimHalf
        else
            dindex(1) = sindex(2)
            dindex(2) = sindex(1) - ndimHalf - 1_IK
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   AIO_ENABLED && XLD_ENABLED && RFP_RDP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim, ndimHalf
        ndim = shape(1)
        ndimHalf = (ndim + 1_IK) / 2_IK
        CHECK_ASSERTION(__LINE__, all(0_IK <= shape), \
        SK_"@getMatIndex(): The condition `all(0_IK <= shape)` must hold. shape, doff = "//getStr(shape))
        CHECK_ASSERTION(__LINE__, shape(1) == shape(2), \
        SK_"@getMatIndex(): The condition `shape(1) == shape(2)` must hold. shape = "//getStr(shape))
        CHECK_ASSERTION(__LINE__, 1_IK <= sindex(2) .and. sindex(2) <= sindex(1), \
        SK_"@getMatIndex(): The condition `1 <= sindex(2) .and. sindex(2) <= sindex(1)` must hold. sindex = "//getStr(sindex))
        CHECK_ASSERTION(__LINE__, all(0_IK < sindex) .and. all(sindex <= shape), \
        SK_"@getMatIndex(): The condition `all(0_IK < sindex) .and. all(sindex <= shape)` must hold. sindex, shape = "//getStr([sindex, shape]))
        if (sindex(2) <= ndimHalf) then
            if (ndim < ndimHalf * 2_IK) then ! odd `ndim`.
                dindex(1) = sindex(1)
            else ! even
                dindex(1) = sindex(1) + 1_IK
            end if
            dindex(2) = sindex(2)
        else
            if (ndim < ndimHalf * 2_IK) then ! odd `ndim`.
                dindex(1) = sindex(2) - ndimHalf
                dindex(2) = sindex(1) - ndimHalf + 1_IK
            else ! even
                dindex(1) = sindex(2) - ndimHalf
                dindex(2) = sindex(1) - ndimHalf
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   AIO_ENABLED && XLD_ENABLED && RDP_RFP_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ndim, ndimHalf, remainder
        CHECK_ASSERTION(__LINE__, all(0_IK <= shape), \
        SK_"@getMatIndex(): The condition `all(0_IK <= shape)` must hold. shape, doff = "//getStr(shape))
        ndim = shape(1)
        ndimHalf = ndim / 2_IK
        remainder = ndim - ndimHalf * 2_IK
        CHECK_ASSERTION(__LINE__, shape(1) == shape(2), \
        SK_"@getMatIndex(): The condition `shape(1) == shape(2)` must hold. shape = "//getStr(shape))
        CHECK_ASSERTION(__LINE__, 1_IK <= sindex(2) .and. sindex(2) <= ndim + 1_IK - remainder, \
        SK_"@getMatIndex(): The condition `1 <= sindex(1) .and. sindex(1) <= ndim + 1 - mod(ndim, 2)` must hold. sindex = "//getStr(sindex))
        CHECK_ASSERTION(__LINE__, 1_IK <= sindex(1) .and. sindex(2) <= ndimHalf + remainder, \
        SK_"@getMatIndex(): The condition `1 <= sindex(2) .and. sindex(2) <= ndimHalf + mod(ndim, 2)` must hold. sindex = "//getStr(sindex))
        if (sindex(2) - sindex(1) < remainder) then
            dindex(1) = sindex(1) + remainder - 1_IK
            dindex(2) = sindex(2)
        else
            dindex(1) = sindex(2) + ndimHalf
            dindex(2) = sindex(1) + ndimHalf + remainder
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
