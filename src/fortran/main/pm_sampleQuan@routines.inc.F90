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
!>  This include file contains procedure implementation of the generic interfaces of [pm_sampleQuan](@ref pm_sampleQuan).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the quantile index.
#if     QD0_ENABLED
#define GET_INDEX(I)I
#elif   QD1_ENABLED
#define GET_INDEX(I):, I
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getQuan_ENABLED && (QD1_ENABLED || QD0_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     ND1_ENABLED
        integer(IK), parameter :: dim = 1
#elif   ND2_ENABLED
        integer(IK) :: idim, ndim, nsam
#else
#error  "Unrecognized interface."
#endif
        real(TKG) :: ecdf(size(sample, dim, IK))
        real(TKG) :: sampleSorted(size(sample, dim, IK))
        CHECK_ASSERTION(__LINE__, all([0._TKG <= prob .and. prob <= 1._TKG]), SK_"@getQuan(): The condition `all([0. <= prob .and. prob <= 1.])` must hold. prob = "//getStr(prob))
        ! Define the optional weight argument.
#if     WNO_ENABLED
        call setECDF(ecdf)
#elif   getQuan_ENABLED && (WTI_ENABLED || WTR_ENABLED)
        if (present(weisum)) then
            call setECDF(ecdf, weight, weisum)
        else
            call setECDF(ecdf, weight, sum(weight))
        end if
#else
#error  "Unrecognized interface."
#endif
        ! Compute the quantiles.
#if     ND1_ENABLED
        sampleSorted = sample
        call setSorted(sampleSorted)
        call setExtrap(method, ecdf, sampleSorted, prob, quan)
#elif   ND2_ENABLED
        CHECK_ASSERTION(__LINE__, dim == 1 .or. dim == 2, SK_"@getQuan(): The condition `dim == 1 .or. dim == 2` must hold. dim = "//getStr(dim))
        ndim = size(sample, 3 - dim, IK)
        nsam = size(sample, dim, IK)
        if (dim == 2) then
            do idim = 1, ndim
                sampleSorted = sample(idim, 1 : nsam)
                call setSorted(sampleSorted)
                call setExtrap(method, ecdf, sampleSorted, prob, quan(GET_INDEX(idim)))
            end do
        else
            do idim = 1, ndim
                sampleSorted = sample(1 : nsam, idim)
                call setSorted(sampleSorted)
                call setExtrap(method, ecdf, sampleSorted, prob, quan(GET_INDEX(idim)))
            end do
        end if
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_INDEX
#undef  WEIGHT