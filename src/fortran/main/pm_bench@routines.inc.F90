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
!>  This file contains procedure implementations of [pm_bench](@ref pm_bench).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Saturday 1:30 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getTiming_ENABLED || setTiming_ENABLED
#if     getTiming_ENABLED
#define TIMING timing
#elif   setTiming_ENABLED
#define TIMING self%timing
#endif
        integer(IK) :: lenTime
        integer(IK) :: counter
        integer(IK) :: miniter_def
        real(RKD)   :: minsec_def
        if (present(miniter)) then
            miniter_def = miniter
        else
            miniter_def = self%miniter
        end if
        if (present(minsec)) then
            minsec_def = minsec
        else
            minsec_def = self%minsec
            !if (present(miniter) .and. miniter_def > 1_IK) then
            !    minsec_def = 0._RKD
            !else
            !    minsec_def = self%minsec
            !end if
        end if
        call self%exec()
        lenTime = 10000_IK
        if (allocated(TIMING%VALUES)) deallocate(TIMING%VALUES)
        allocate(TIMING%VALUES(lenTime))
        counter = 1_IK
        TIMING%overhead = 0._RKD
        TIMING%mean = 0._RKD
        loopRepeatTiming: do
            self%timer%start = self%timer%time()
            call self%exec()
            TIMING%VALUES(counter) = max(self%timer%time(since = self%timer%start), self%timer%resol)
            TIMING%mean = TIMING%mean + TIMING%VALUES(counter)
            self%timer%start = self%timer%time()
            call self%overhead()
            TIMING%overhead = TIMING%overhead + self%timer%time(since = self%timer%start)
            if (minsec_def <= TIMING%mean .and. miniter_def <= counter) exit loopRepeatTiming
            if (counter == lenTime) then
                lenTime = 2_IK * lenTime
                call setResized(TIMING%VALUES, size = lenTime)
            end if
            counter = counter + 1_IK
        end do loopRepeatTiming
        if (size(TIMING%VALUES, kind = IK) /= counter) call setResized(TIMING%VALUES, size = counter)
        TIMING%overhead = TIMING%overhead / counter
        TIMING%mean = max(self%timer%resol, TIMING%mean / counter - TIMING%overhead)
        TIMING%std = 0._RKC
        if (1_IK < size(TIMING%VALUES, 1, IK)) TIMING%std = max(self%timer%resol, sqrt(getVar(TIMING%VALUES)))
        TIMING%min = max(self%timer%resol, minval(TIMING%VALUES) - TIMING%overhead)
        TIMING%max = max(self%timer%resol, maxval(TIMING%VALUES) - TIMING%overhead)
#else
#error  "Unrecognized interface."
#endif
#undef  TIMING