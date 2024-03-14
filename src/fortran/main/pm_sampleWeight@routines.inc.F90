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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_sampleWeight](@ref pm_sampleWeight).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Saturday 2:33 AM, August 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK_ENABLED
#define TYPE_OF_SKIP integer(IKC)
#elif   RK_ENABLED || IK_RK_ENABLED
#define TYPE_OF_SKIP real(RKC)
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%
#if     getReweight_ENABLED
        !%%%%%%%%%%%%%%%%%%

        reweight = weight
        call setReweight(reweight, skip)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setReweight_ENABLED && (IK_ENABLED || IK_RK_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_SKIP :: current
        integer(IK) :: isam
        current = 0
        CHECK_ASSERTION(__LINE__, 0 < skip, SK_": The condition `0 < skip` must hold. skip = "//getStr(skip))
        loopOverWeight: do isam = 1, size(weight, 1, IK)
            if (0_IKC < weight(isam)) then
                current = current + weight(isam)
                weight(isam) = 0_IKC
                loopOverElement: do
                    if (current < skip) exit loopOverElement
                    weight(isam) = weight(isam) + 1_IKC
                    current = current - skip
                end do loopOverElement
            else
                weight(isam) = 0_IKC
            end if
        end do loopOverWeight

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setReweight_ENABLED && RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_SKIP :: current
        integer(IK) :: isam
        current = 0
        CHECK_ASSERTION(__LINE__, 0 < skip, SK_": The condition `0 < skip` must hold. skip = "//getStr(skip))
        loopOverWeight: do isam = 1, size(weight, 1, IK)
            if (0._RKC < weight(isam)) then
                current = current + weight(isam)
                weight(isam) = 0._RKC
                loopOverElement: do
                    if (current < skip) exit loopOverElement
                    weight(isam) = weight(isam) + 1._RKC
                    current = current - skip
                end do loopOverElement
            else
                weight(isam) = 0._RKC
            end if
        end do loopOverWeight
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_SKIP