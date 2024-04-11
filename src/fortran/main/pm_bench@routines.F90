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

submodule (pm_bench) routines ! LCOV_EXCL_LINE

    use pm_sampleVar, only: getVar
    use pm_arrayResize, only: setResized
    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure benchBase_typer
        use pm_timer, only: timer_type
        benchBase%name = trim(adjustl(name))
        if (present(minsec)) benchBase%minsec = minsec
        if (present(miniter)) benchBase%miniter = miniter
        if (present(timer)) then
            allocate(benchBase%timer, source = timer)
        else
            allocate(benchBase%timer, source = timer_type())
        end if
        !benchBase%timing%overhead = benchBase%timer%resol ! minimum overhead
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure bench_typer
        bench%benchBase_type = benchBase_type(name, minsec = minsec, miniter = miniter, timer = timer)
        if (present(overhead)) then
            bench%overhead => overhead
        else
            bench%overhead => doNothing
        end if
        bench%exec => exec
        if (present(minsec)) bench%minsec = minsec
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure getTimingMethod
        use pm_kind, only: RKC => RKD
#define getTiming_ENABLED 1
#include "pm_bench@routines.inc.F90"
#undef  getTiming_ENABLED
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure setTimingMethod
        use pm_kind, only: RKC => RKD
#define setTiming_ENABLED 1
#include "pm_bench@routines.inc.F90"
#undef  setTiming_ENABLED
    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure benchMulti_typer

        use pm_kind, only: IK, LK
        use pm_arrayRange, only: setRange
        use pm_arrayShuffle, only: setShuffled
        use pm_arrayReverse, only: setReversed

        character(*, SK), parameter :: namesep = SK_" vs. "
        integer(IK) , allocatable :: index(:)
        integer(IK) :: repeat_def
        logical(LK) :: sorted_def
        integer(IK) :: i, lenCase
        integer(IK) :: lenCaseTot

        self%ncase = size(case, 1, IK)
        if (0_IK < self%ncase) then
            self%name = case(1)%name
            do i = 2, self%ncase
                self%name = self%name//namesep//case(i)%name
            end do
        else
            self%name = ""
        end if

        repeat_def = 2_IK
        sorted_def = .false._LK
        if (present(sorted)) sorted_def = sorted
        if (present(repeat)) repeat_def = repeat
        lenCase = size(case, kind = IK)
        allocate(index(size(case, kind = IK) * repeat_def))
        call setRange(index(1:lenCase), 1_IK)
        do i = 1_IK, repeat_def - 1_IK
            call setReversed(index((i-1) * lenCase + 1 : i * lenCase), index(i * lenCase + 1 : (i + 1) * lenCase))
        end do
        if (.not. sorted_def) call setShuffled(index)
        lenCaseTot = lenCase * repeat_def
        allocate(self%case(lenCaseTot))
        do i = 1_IK, lenCaseTot
            self%case(i) = case(index(i))
            call self%case(i)%setTiming()
        end do
        deallocate(index)

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module procedure showsum_

        use pm_kind, only: RKD
        use pm_kind, only: IK, LK
        use pm_val2str, only: getStr
        use iso_fortran_env, only: output_unit
        use pm_arraySort, only: setSorted

        logical(LK) :: tabular_def
        integer(IK) :: unit_def, lenNameMax, i
        integer(IK) :: rank(size(self%case))
        real(RKD)   :: meanRunTime(size(self%case))
        character(:, SK), allocatable :: format

        unit_def = output_unit
        if (present(unit)) unit_def = unit

        if (present(tabular)) then
            tabular_def = tabular
        elseif (unit_def == output_unit) then
            tabular_def = .true._LK ! LCOV_EXCL_LINE
        else
            tabular_def = .false._LK
        end if

        ! Fetch the mean runtimes and determine the maximum benchmark name length.

        lenNameMax = 0_IK
        do i = 1_IK, size(self%case, kind = IK)
            meanRunTime(i) = self%case(i)%timing%mean
            if (lenNameMax < len(self%case(i)%name, IK)) lenNameMax = len(self%case(i)%name, IK)
        end do

        ! Get the sorted indices of the mean runtimes.

        call setSorted(meanRunTime, rank)
        if (tabular_def) then
            format = SK_"(A"//getStr(lenNameMax)//SK_",' | ', *(g15.9,:,' | '))"
            write(unit_def, *)
            write(unit_def, format) "benchmark", "avg runtime [s]", "std runtime [s]"
            write(unit_def, format) repeat(SK_"-", lenNameMax), repeat(SK_"-", 15), repeat(SK_"-", 15)
            do i = 1_IK, size(meanRunTime, 1, IK)
                write(unit_def, format) self%case(rank(i))%name, & ! LCOV_EXCL_LINE
                                        self%case(rank(i))%timing%mean, & ! LCOV_EXCL_LINE
                                        self%case(rank(i))%timing%std
            end do
            write(unit_def, *)
        else
            format = SK_"(*(g0.9,:,','))"
            write(unit_def, format) "benchmark", "avg runtime [s]", "std runtime [s]"
            do i = 1_IK, size(meanRunTime, 1, IK)
                write(unit_def, format) self%case(rank(i))%name, & ! LCOV_EXCL_LINE
                                        self%case(rank(i))%timing%mean, & ! LCOV_EXCL_LINE
                                        self%case(rank(i))%timing%std
            end do
        end if

    end procedure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule routines