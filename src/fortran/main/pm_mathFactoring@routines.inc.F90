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
!>  This include file contains implementations of the procedures in module [pm_mathFactoring](@ref pm_mathFactoring).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getFactoring_IK_ENABLED
        use pm_arrayResize, only: setResized
        integer(IK)     :: count, csize
        integer(IKG)    :: halfn, divisor
        CHECK_ASSERTION(__LINE__, 1_IKG < posint, \
        SK_"@getFactoring(): The condition `1 < posint` must hold for corresponding input arguments. posint = "//getStr(posint))
        csize = 15_IK
        allocate(Factoring(csize))
        count = 0_IK
        do ! first remove all factors of 2.
            halfn = posint / 2_IKG
            if (halfn * 2_IKG /= posint) exit
            count = count + 1_IK
            if (csize < count) then
                csize = csize * 2_IK
                call setResized(Factoring, csize)
            end if
            Factoring(count) = 2_IKG
            posint = halfn
        end do
        ! Find the odd factors.
        divisor = 3_IK
        do  !   3, 5, 7, .... will be tried. 
            !   \todo This algorithm can be improved.
            if (divisor > posint) exit ! If a factor is too large, we are done.
            do  ! Try the current factor repeatedly, until all is taken out.
                if (mod(posint, divisor) /= 0_IKG .or. posint == 1_IKG)  exit
                count = count + 1_IK
                if (csize < count) then
                    csize = csize * 2_IK
                    call setResized(Factoring, csize)
                end if
                Factoring(count) = divisor
                posint = posint / divisor ! Remove the current factor from `posint`.
            end do
            divisor = divisor + 2_IKG ! Move to next odd number.
        end do
        Factoring = Factoring(1:count)
#else
#error  "Unrecognized interface."
#endif