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
!>  This file contains procedure implementations of [pm_arrayCopy](@ref pm_arrayCopy).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     D0_ENABLED && SK_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE(x) len(x, kind = IK)
#elif   D1_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE(x) size(x, kind = IK)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#if     setCopyIndexed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i
        CHECK_ASSERTION(__LINE__, size(indexF, 1, IK) == size(indexT, 1, IK), \
        SK_"@setCopyIndexed(): The condition `size(indexF) == size(indexT)` must hold. size(indexF), size(indexT) = "\
        //getStr([size(indexF, 1, IK) == size(indexT, 1, IK)])) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK < indexF) .and. all(indexF <= GET_SIZE(From)), \
        SK_"@setCopyIndexed(): The condition `all(1 < indexF) .and. all(indexF <= size(From))` must hold. size(From), indexF = "\
        //getStr([GET_SIZE(From), indexF])) ! fpp
        CHECK_ASSERTION(__LINE__, all(0_IK < indexT) .and. all(indexT <= GET_SIZE(To  )), \
        SK_"@setCopyIndexed(): The condition `all(1 < indexT) .and. all(indexT <= size(To  ))` must hold. size(To  ), indexT = "\
        //getStr([GET_SIZE(To  ), indexT])) ! fpp
        do i = 1_IK, size(indexF, 1, IK)
            To(GET_INDEX(indexT(i))) = From(GET_INDEX(indexF(i)))
        end do

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setCopyStrided_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: ifrom, ito
#if     CHECK_ENABLED
        if (incf /= 0_IK .and. inct /= 0_IK) CHECK_ASSERTION(__LINE__, (GET_SIZE(From) - 1_IK) / abs(incf) == (GET_SIZE(To) - 1_IK) / abs(inct), \
        SK_"@setCopyStrided(): The condition `(size(From)-1)/abs(incf) == (size(To)-1)/abs(inct)` must hold. size(From), size(To), incf, inct = "\
        //getStr([GET_SIZE(From), GET_SIZE(To), incf, inct]))
#endif
        if (incf > 0_IK) then
            ito = 1_IK
            if (inct < 0_IK) ito = GET_SIZE(To)
            do ifrom = 1_IK, GET_SIZE(From), incf
                To(GET_INDEX(ito)) = From(GET_INDEX(ifrom))
                ito = ito + inct
            end do
        elseif (incf < 0_IK) then
            ito = 1_IK
            if (inct < 0_IK) ito = GET_SIZE(To)
            do ifrom = GET_SIZE(From), 1_IK, incf
                To(GET_INDEX(ito)) = From(GET_INDEX(ifrom))
                ito = ito + inct
            end do
        elseif (inct > 0_IK) then
            do concurrent(ito = 1_IK : GET_SIZE(To) : inct)
                To(GET_INDEX(ito)) = From(GET_INDEX(1_IK))
            end do
        elseif (inct < 0_IK) then
            do concurrent(ito = GET_SIZE(To) : 1_IK : inct)
                To(GET_INDEX(ito)) = From(GET_INDEX(1_IK))
            end do
        else
            error stop "The condition `incf /= 0_IK .or. inct /= 0_IK` must hold." ! LCOV_EXCL_LINE
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef GET_INDEX
#undef GET_SIZE