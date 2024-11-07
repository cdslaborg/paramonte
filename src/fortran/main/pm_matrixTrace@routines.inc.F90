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
!>  This file contains procedure implementations of [pm_matrixTrace](@ref pm_matrixTrace).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the increment operations for various definitions of trace.
#if     getMatTrace_ENABLED
#define TRACE trace
#define OPERATION +
#define GET_LOG(X) X
#elif   getMatMulTrace_ENABLED
#define TRACE mulTrace
#define OPERATION *
#define GET_LOG(X) X
#elif   getMatMulTraceLog_ENABLED
#define TRACE logMulTrace
#define OPERATION +
#if     IK_ENABLED
#define GET_LOG(X) log(real(X, RKD))
#elif   CK_ENABLED || RK_ENABLED
#define GET_LOG(X) log(X)
#else
#error  "Unrecognized interface."
#endif
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     (getMatTrace_ENABLED || getMatMulTrace_ENABLED || getMatMulTraceLog_ENABLED) && (DEF_ENABLED || RDP_ENABLED) && XXX_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim
        CHECK_ASSERTION(__LINE__, all(0_IK < shape(mat, IK)), SK_"@getMatTrace(): The condition `all(0 < shape(mat, IK))` must hold. shape(mat) = "//getStr(shape(mat, IK)))
        CHECK_ASSERTION(__LINE__, size(mat, 1, IK) == size(mat, 2, IK), SK_"@getMatTrace(): The condition `size(mat, 1) == size(mat, 2)` must hold. shape(mat) = "//getStr(shape(mat, IK)))
        TRACE = GET_LOG(mat(1, 1))
        do idim = 2, size(mat, 1, IK)
            TRACE = TRACE OPERATION GET_LOG(mat(idim, idim))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getMatTrace_ENABLED || getMatMulTrace_ENABLED || getMatMulTraceLog_ENABLED) && RFP_ENABLED && (UXD_ENABLED || XLD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, nrow, ncol
        nrow = size(mat, 1, IK)
        ncol = size(mat, 2, IK)
        CHECK_ASSERTION(__LINE__, all(0_IK < shape(mat, IK)), SK_"@getMatTrace(): The condition `all(0 < shape(mat, IK))` must hold. shape(mat) = "//getStr(shape(mat, IK)))
        if (ncol < nrow) then ! regular RFP packing, no transposition.
            if (nrow < 2 * ncol) then ! odd order.
                CHECK_ASSERTION(__LINE__, nrow == 2 * ncol - 1, SK_"@getMatTrace(): The condition `size(mat, 1) == 2 * size(mat, 2) - 1)` must hold. shape(mat) = "//getStr(shape(mat, IK)))
#if             UXD_ENABLED
#define         OFFSET_PLUS(X) X +
                TRACE = GET_LOG(mat(ncol, 1))
                do idim = 1, ncol - 1
                    TRACE = TRACE OPERATION GET_LOG(mat(ncol + idim, idim))
                    TRACE = TRACE OPERATION GET_LOG(mat(ncol + idim, idim + 1))
                end do
#elif           XLD_ENABLED
#define         OFFSET_PLUS(X)
                TRACE = GET_LOG(mat(1, 1))
                do idim = 2, ncol
                    TRACE = TRACE OPERATION GET_LOG(mat(idim - 1, idim))
                    TRACE = TRACE OPERATION GET_LOG(mat(idim, idim))
                end do
#endif
            else ! even order.
                CHECK_ASSERTION(__LINE__, nrow == 2 * ncol + 1, SK_"@getMatTrace(): The condition `size(mat, 1) == 2 * size(mat, 2) + 1)` must hold. shape(mat) = "//getStr(shape(mat, IK)))
                idim = 1
                TRACE = GET_LOG(mat(OFFSET_PLUS(ncol) idim, 1))
                TRACE = TRACE OPERATION GET_LOG(mat(OFFSET_PLUS(ncol) idim + 1, 1))
                do idim = 2, ncol
                    TRACE = TRACE OPERATION GET_LOG(mat(OFFSET_PLUS(ncol) idim      , idim))
                    TRACE = TRACE OPERATION GET_LOG(mat(OFFSET_PLUS(ncol) idim + 1  , idim))
                end do
            end if
        elseif (nrow < ncol) then ! transposed RFP packing.
            if (ncol < 2 * nrow) then ! odd order.
                CHECK_ASSERTION(__LINE__, ncol == 2 * nrow - 1, SK_"@getMatTrace(): The condition `size(mat, 2) == 2 * size(mat, 1) - 1)` must hold. shape(mat) = "//getStr(shape(mat, IK)))
#if             UXD_ENABLED
                TRACE = GET_LOG(mat(1, nrow))
                do idim = 1, nrow - 1
                    TRACE = TRACE OPERATION GET_LOG(mat(idim    , nrow + idim))
                    TRACE = TRACE OPERATION GET_LOG(mat(idim + 1, nrow + idim))
                end do
#elif           XLD_ENABLED
                TRACE = GET_LOG(mat(1, 1))
                do idim = 2, nrow
                    TRACE = TRACE OPERATION GET_LOG(mat(idim, idim - 1))
                    TRACE = TRACE OPERATION GET_LOG(mat(idim, idim))
                end do
#endif
            else ! even order.
                CHECK_ASSERTION(__LINE__, ncol == 2 * nrow + 1, SK_"@getMatTrace(): The condition `size(mat, 2) == 2 * size(mat, 1) + 1)` must hold. shape(mat) = "//getStr(shape(mat, IK)))
                idim = 1
                TRACE = GET_LOG(mat(1, OFFSET_PLUS(nrow) idim))
                TRACE = TRACE OPERATION GET_LOG(mat(1, OFFSET_PLUS(nrow) idim + 1))
                do idim = 2, nrow
                    TRACE = TRACE OPERATION GET_LOG(mat(idim, OFFSET_PLUS(nrow) idim    ))
                    TRACE = TRACE OPERATION GET_LOG(mat(idim, OFFSET_PLUS(nrow) idim + 1))
                end do
            end if
        else ! nrow == ncol == 1
            CHECK_ASSERTION(__LINE__, all(1_IK == shape(mat, IK)), SK_"@getMatTrace(): The input RFP `mat` must conform to an upper-triangular subset in RDP packing. shape(mat) = "//getStr(shape(mat, IK)))
            TRACE = GET_LOG(mat(1, 1))
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   (getMatTrace_ENABLED || getMatMulTrace_ENABLED || getMatMulTraceLog_ENABLED) && LFP_ENABLED && (UXD_ENABLED || XLD_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: idim, index, ndim, nell
        nell = size(mat, 1, IK)
        CHECK_ASSERTION(__LINE__, 0_IK < nell, SK_"@getMatTrace(): The condition `0 < size(mat)` must hold. size(mat) = "//getStr(nell))
        !ndim = (getSqrt(nell * 8 + 1) - 1) / 2
        ndim = (getSqrt(nell * 8 + 1) - 1) / 2
        CHECK_ASSERTION(__LINE__, nell == ndim * (ndim + 1) / 2, SK_"@getMatTrace(): The input LFP `mat` must conform to an upper-triangular subset in RDP packing. ndim, shape(mat) = "//getStr([ndim, shape(mat, IK)]))
        TRACE = GET_LOG(mat(1))
        if (1_IK < nell) then
#if         UXD_ENABLED
            index = 1
            do idim = 2, ndim
                index = index + idim
                TRACE = TRACE OPERATION GET_LOG(mat(index))
            end do
#elif       XLD_ENABLED
            index = 0
            do idim = 0, ndim - 2
                index = index + ndim - idim
                TRACE = TRACE OPERATION GET_LOG(mat(index))
            end do
#endif
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  OFFSET_PLUS
#undef  OPERATION
#undef  GET_LOG
#undef  TRACE